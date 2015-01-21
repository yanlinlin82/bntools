#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "base_map.h"
#include "nick_map.h"

#define DEF_ENZ_NAME "BspQI"
#define DEF_REC_SEQ "GCTCTTCN^"
#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "bed"

#define MAX_ENZYME_NAME_SIZE 16
#define MAX_REC_SEQ_SIZE 32
#define MAX_CHROM_NAME_SIZE 64

static int verbose = 0;

static char enzyme_name[MAX_ENZYME_NAME_SIZE] = DEF_ENZ_NAME;
static char rec_seq[MAX_REC_SEQ_SIZE] = DEF_REC_SEQ;
static char output_file[PATH_MAX] = DEF_OUTPUT;
static int output_cmap = 0;
static int transform_to_number = 0;

static char rec_bases[MAX_REC_SEQ_SIZE] = { };
static int rec_seq_size = 0;
static int nick_offset = -1;
static int palindrome = 1;

static char buf[MAX_REC_SEQ_SIZE * 2];
static int pos = 0;

struct nick_map map = { };

static int prepare_rec_seq(void)
{
	int i;
	for (i = 0, rec_seq_size = 0; rec_seq[i]; ++i) {
		if (rec_seq[i] == '^') {
			if (nick_offset >= 0) {
				fprintf(stderr, "Error: Invalid recognition sequence '%s'\n", rec_seq);
				return 1;
			}
			nick_offset = i;
		} else {
			char c = BASE_MAP[(int)rec_seq[i]];
			if (c == 0) {
				fprintf(stderr, "Error: Invalid character '%c' in recognition sequence '%s'\n", rec_seq[i], rec_seq);
				return 1;
			}
			if (rec_seq_size >= sizeof(rec_bases)) {
				fprintf(stderr, "Error: recognition sequence is too long\n");
				return 1;
			}
			rec_bases[rec_seq_size++] = c;
		}
	}
	if (nick_offset < 0) {
		nick_offset = 0;
	}
	for (i = 0; i < rec_seq_size / 2; ++i) {
		if (rec_bases[i] != (~rec_bases[rec_seq_size - i - 1] & BASE_N)) {
			palindrome = 0;
			break;
		}
	}
	return 0;
}

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools nick [options] <x.fa>\n"
			"\n"
			"Options:\n"
			"   <x.fa>          input FASTA sequence to generate restriction map\n"
			"   -v              show verbose messages\n"
			"   -o FILE         output file ["DEF_OUTPUT"]\n"
			"   -f {bed|cmap}   output format ["DEF_FORMAT"]\n"
			"   -e STR          restriction enzyme name ["DEF_ENZ_NAME"]\n"
			"   -r STR          recognition sequence ["DEF_REC_SEQ"]\n"
			"   -n              transform chromosome names into numbers\n"
			"\n");
}

static inline int is_whitespace(int c)
{
	return (c == ' ' || c == '\t' || c == '\r' || c == '\n');
}

static int seq_match(const char *ref, const char *query, size_t len, int revcomp)
{
	size_t i;
	for (i = 0; i < len; ++i) {
		int c;
		switch (ref[i]) {
		case BASE_A: c = (revcomp ? BASE_T : BASE_A); break;
		case BASE_C: c = (revcomp ? BASE_G : BASE_C); break;
		case BASE_G: c = (revcomp ? BASE_C : BASE_G); break;
		case BASE_T: c = (revcomp ? BASE_A : BASE_T); break;
		default: return 0;
		}
		if ((c & query[revcomp ? (len - i - 1) : i]) != c) {
			return 0;
		}
	}
	return 1;
}

static int process_line(struct nick_site_list *list, const char *line, const char *chrom, int base_count)
{
	const char *p;
	int revcomp;
	for (p = line; *p; ++p) {
		if (is_whitespace(*p)) {
			continue;
		}
		if (pos >= sizeof(buf)) {
			memcpy(buf, buf + sizeof(buf) - rec_seq_size + 1, rec_seq_size - 1);
			pos = rec_seq_size - 1;
		}
		++base_count;
		buf[pos++] = BASE_MAP[(int)*p];
		if (pos < rec_seq_size) {
			continue;
		}
		for (revcomp = 0; revcomp <= (palindrome ? 0 : 1); ++revcomp) {
			if (seq_match(buf + pos - rec_seq_size, rec_bases, rec_seq_size, revcomp)) {
				int site_pos = base_count - (revcomp ? nick_offset : (rec_seq_size - nick_offset));
				if (nick_map_add_site(list, site_pos, revcomp)) {
					return -ENOMEM;
				}
			}
		}
	}
	return base_count;
}

static void write_cmap_header(gzFile fout, size_t seq_total_number)
{
	gzprintf(fout, "# CMAP File Version:  0.1\n");
	gzprintf(fout, "# Label Channels:  1\n");
	gzprintf(fout, "# Nickase Recognition Site 1:  %s/%s\n", enzyme_name, rec_seq);
	gzprintf(fout, "# Number of Consensus Nanomaps:    %zd\n", seq_total_number);
	gzprintf(fout, "#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence\n");
	gzprintf(fout, "#f int	float	int	int	int	float	float	int	int\n");
}

static void output_cmap_line(gzFile fout, const char *cmap_id, int contig_length,
		size_t num_sites, size_t site_id, int label_channel, int position,
		int stddev, int coverage, int occurance)
{
	gzprintf(fout, "%s\t%d\t%zd\t%zd\t%d\t%d\t%d\t%d\t%d\n",
			cmap_id, contig_length, num_sites, site_id,
			label_channel, position, stddev, coverage, occurance);
}

static void print_nick_site_list(gzFile fout, const struct nick_site_list *p)
{
	size_t j, k, count;
	if (output_cmap) {
		for (j = 0, count = 0; j < p->size; ++j) {
			if (j > 0 && p->data[j - 1].pos == p->data[j].pos) {
				continue;
			}
			++count;
		}
		for (j = 0, k = 0; j < p->size; ++j) {
			if (j > 0 && p->data[j - 1].pos == p->data[j].pos) {
				continue;
			}
			output_cmap_line(fout, p->chrom_name, p->chrom_size, count,
					++k, 1, p->data[j].pos, 0, 0, 0);
		}
		if (count > 0) {
			output_cmap_line(fout, p->chrom_name, p->chrom_size, count,
					count + 1, 0, p->chrom_size, 0, 1, 1);
		}
	} else {
		for (j = 0, k = 0; j < p->size; ++j) {
			gzprintf(fout, "%s\t%d\t%d\t%s/%s\t0\t%c\n",
					p->chrom_name, p->data[j].pos, p->data[j].pos + 1,
					enzyme_name, rec_seq, "+-"[p->data[j].strand]);
		}
	}
}

static void print_results(gzFile fout)
{
	size_t i;
	if (output_cmap) {
		write_cmap_header(fout, map.size);
	}
	for (i = 0; i < map.size; ++i) {
		print_nick_site_list(fout, &map.data[i]);
	}
}

static int process(gzFile fin, gzFile fout)
{
	struct nick_site_list *list = NULL;
	char chrom[MAX_CHROM_NAME_SIZE] = "";
	int base_count = 0;
	while (!gzeof(fin)) {
		char buf[256];
		if (!gzgets(fin, buf, sizeof(buf))) break;
		if (buf[0] == '>') {
			if (list) {
				if (!output_cmap) {
					print_nick_site_list(fout, list);
				}
				list->chrom_size = base_count;
			}
			if (transform_to_number) {
				static int number = 0;
				snprintf(chrom, sizeof(chrom), "%d", ++number);
			} else {
				char *p = buf + 1;
				while (*p && !isspace(*p)) ++p;
				*p = '\0';
				snprintf(chrom, sizeof(chrom), buf + 1);
			}
			if (verbose > 0) {
				fprintf(stderr, "Loading sequence '%s' ... ", chrom);
			}
			base_count = 0;
			list = nick_map_add_chrom(&map, chrom);
			if (!list) {
				return -ENOMEM;
			}
		} else {
			int n = process_line(list, buf, chrom, base_count);
			if (n < 0) {
				return n;
			}
			base_count = n;
		}
	}
	if (list) {
		list->chrom_size = base_count;
		if (!output_cmap) {
			print_nick_site_list(fout, list);
		}
	}
	if (output_cmap) {
		print_results(fout);
	}
	return 0;
}

static int nick(const char *in)
{
	gzFile fin = NULL;
	gzFile fout = NULL;
	int c;

	if (strcmp(in, "-") == 0 || strcmp(in, "stdin") == 0) {
		fin = gzdopen(0, "r"); /* stdin */
	} else {
		fin = gzopen(in, "r");
	}
	if (!fin) {
		fprintf(stderr, "Error: Can not open FASTA file '%s'\n", in);
		return 1;
	}

	if (strcmp(output_file, "-") == 0 || strcmp(output_file, "stdout") == 0) {
		fout = gzdopen(1, "wT"); /* stdout, without compression */
	} else {
		size_t len = strlen(output_file);
		if (len > 3 && strcmp(output_file + len - 3, ".gz") == 0) {
			fout = gzopen(output_file, "wx"); /* 'x' is for checking existance */
		} else {
			fout = gzopen(output_file, "wxT"); /* without compression */
		}
	}
	if (!fout) {
		if (errno == EEXIST) {
			fprintf(stderr, "Error: Output file '%s' has already existed!\n", output_file);
		} else {
			fprintf(stderr, "Error: Can not open output file '%s'\n", output_file);
		}
		gzclose(fin);
		return 1;
	}

	c = gzgetc(fin);
	if (c != '>') {
		fprintf(stderr, "Error: File '%s' is not in FASTA format\n", in);
		gzclose(fout);
		gzclose(fin);
		return 1;
	}
	gzungetc(c, fin);

	if (process(fin, fout)) {
		return 1;
	}

	gzclose(fout);
	gzclose(fin);
	nick_map_free(&map);
	return 0;
}

int nick_main(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "vo:f:e:r:n")) != -1) {
		switch (c) {
		case 'v':
			++verbose;
			break;
		case 'o':
			snprintf(output_file, sizeof(output_file), "%s", optarg);
			break;
		case 'f':
			if (strcmp(optarg, "bed") == 0) {
				output_cmap = 0;
			} else if (strcmp(optarg, "cmap") == 0) {
				output_cmap = 1;
			} else {
				fprintf(stderr, "Error: Unknown output format '%s'!\n", optarg);
				return 1;
			}
			break;
		case 'e':
			snprintf(enzyme_name, sizeof(enzyme_name), "%s", optarg);
			break;
		case 'r':
			snprintf(rec_seq, sizeof(rec_seq), "%s", optarg);
			break;
		case 'n':
			transform_to_number = 1;
			break;
		default:
			return 1;
		}
	}
	if (argc != optind + 1) {
		print_usage();
		return 1;
	}

	init_base_map();

	if (prepare_rec_seq() != 0) {
		return 1;
	}
	return nick(argv[optind]);
}
