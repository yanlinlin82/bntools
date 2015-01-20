#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "nick.h"

#define NAME "nick"
#define DEF_ENZ_NAME "BspQI"
#define DEF_REC_SEQ "GCTCTTCN^"
#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "bed"

#define MAX_ENZYME_NAME_SIZE 16
#define MAX_REC_SEQ_SIZE 32
#define MAX_CHROM_NAME_SIZE 64

#define BASE_A 1
#define BASE_C 2
#define BASE_G 4
#define BASE_T 8
#define BASE_N (BASE_A | BASE_C | BASE_G | BASE_T)

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

static char chrom[MAX_CHROM_NAME_SIZE] = "";
static int offset = 0;

static char buf[MAX_REC_SEQ_SIZE * 2];
static int pos = 0;

struct nick_site {
	int pos;
	int revcomp;
};

struct nick_site_list {
	size_t capacity;
	size_t size;
	struct nick_site *data;
};

struct nick_site_list results = { };

#define MIN_INCREMENT 16
#define MAX_INCREMENT 1024

static size_t new_capacity(size_t capacity)
{
	if (capacity < MIN_INCREMENT) {
		return capacity + MIN_INCREMENT;
	} else if (capacity > MAX_INCREMENT) {
		return capacity + MAX_INCREMENT;
	} else {
		return capacity * 2;
	}
}

static void add_result(int pos, int revcomp)
{
	size_t i;

	if (results.size + 1 >= results.capacity) {
		size_t capacity = new_capacity(results.capacity);
		struct nick_site *p = malloc(sizeof(struct nick_site) * capacity);
		if (!p) {
			fprintf(stderr, "Error: out of memory");
			exit(1);
		}
		if (results.size > 0) {
			memcpy(p, results.data, sizeof(struct nick_site) * results.size);
		}
		free(results.data);
		results.data = p;
		results.capacity = capacity;
	}
	assert(results.size + 1 < results.capacity);

	for (i = results.size; i > 0; --i) {
		if (results.data[i - 1].pos > pos) {
			memcpy(&results.data[i], &results.data[i - 1], sizeof(struct nick_site));
		} else {
			break;
		}
	}
	results.data[i].pos = pos;
	results.data[i].revcomp = revcomp;
	results.size++;
}

static void free_results(void)
{
	free(results.data);
	results.data = NULL;
	results.size = 0;
	results.capacity = 0;
}

static char BASE_MAP[256] = { };
static char BASE_REV_MAP[16] = { };

static void init_base_map(void)
{
	BASE_MAP['A'] = BASE_MAP['a'] = BASE_A;
	BASE_MAP['C'] = BASE_MAP['c'] = BASE_C;
	BASE_MAP['G'] = BASE_MAP['g'] = BASE_G;
	BASE_MAP['T'] = BASE_MAP['t'] = BASE_T;
	BASE_MAP['U'] = BASE_MAP['u'] = BASE_T;
	BASE_MAP['M'] = BASE_MAP['m'] = BASE_A | BASE_C; /* aMino */
	BASE_MAP['K'] = BASE_MAP['k'] = BASE_G | BASE_T; /* Keto */
	BASE_MAP['R'] = BASE_MAP['r'] = BASE_A | BASE_G; /* puRine */
	BASE_MAP['Y'] = BASE_MAP['y'] = BASE_C | BASE_T; /* pYrimidine */
	BASE_MAP['S'] = BASE_MAP['s'] = BASE_C | BASE_G; /* strong */
	BASE_MAP['W'] = BASE_MAP['w'] = BASE_A | BASE_T; /* weak */
	BASE_MAP['B'] = BASE_MAP['b'] = BASE_C | BASE_G | BASE_T; /* not 'A' */
	BASE_MAP['D'] = BASE_MAP['d'] = BASE_A | BASE_G | BASE_T; /* not 'C' */
	BASE_MAP['H'] = BASE_MAP['h'] = BASE_A | BASE_C | BASE_T; /* not 'G' */
	BASE_MAP['V'] = BASE_MAP['v'] = BASE_A | BASE_C | BASE_G; /* not 'T/U' */
	BASE_MAP['N'] = BASE_MAP['n'] = BASE_N;
	BASE_MAP['X'] = BASE_MAP['x'] = BASE_N;

	BASE_REV_MAP[BASE_A] = 'A';
	BASE_REV_MAP[BASE_C] = 'C';
	BASE_REV_MAP[BASE_G] = 'G';
	BASE_REV_MAP[BASE_T] = 'T';
	BASE_REV_MAP[BASE_A | BASE_C] = 'M';
	BASE_REV_MAP[BASE_G | BASE_T] = 'K';
	BASE_REV_MAP[BASE_A | BASE_G] = 'R';
	BASE_REV_MAP[BASE_C | BASE_T] = 'Y';
	BASE_REV_MAP[BASE_C | BASE_G] = 'S';
	BASE_REV_MAP[BASE_A | BASE_T] = 'W';
	BASE_REV_MAP[BASE_C | BASE_G | BASE_T] = 'B';
	BASE_REV_MAP[BASE_A | BASE_G | BASE_T] = 'D';
	BASE_REV_MAP[BASE_A | BASE_C | BASE_T] = 'H';
	BASE_REV_MAP[BASE_A | BASE_C | BASE_G] = 'V';
	BASE_REV_MAP[BASE_N] = 'N';
}

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
			"Usage: bntools "NAME" [options] <x.fa>\n"
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

static void process_line(gzFile fout, const char *line, const char *chrom)
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
		++offset;
		buf[pos++] = BASE_MAP[(int)*p];
		if (pos < rec_seq_size) {
			continue;
		}
		for (revcomp = 0; revcomp <= (palindrome ? 0 : 1); ++revcomp) {
			if (seq_match(buf + pos - rec_seq_size, rec_bases, rec_seq_size, revcomp)) {
				int pos = offset - (revcomp ? nick_offset : (rec_seq_size - nick_offset));
				add_result(pos, revcomp);
			}
		}
	}
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

static void print_results(gzFile fout, const char *chrom, int chrom_size)
{
	size_t i;
	for (i = 0; i < results.size; ++i) {
		const struct nick_site *p = &results.data[i];
		if (output_cmap) {
			output_cmap_line(fout, chrom, chrom_size, results.size,
					i + 1, 1, p->pos, 0, 0, 0);
		} else {
			gzprintf(fout, "%s\t%d\t%d\t%s/%s\t0\t%c\n",
					chrom, p->pos, p->pos + 1, enzyme_name, rec_seq,
					(p->revcomp ? '-' : '+'));
		}
	}
	if (output_cmap && results.size > 0) {
		output_cmap_line(fout, chrom, chrom_size, results.size,
				results.size + 1, 0, chrom_size, 0, 1, 1);
	}
	results.size = 0;
}

static size_t get_seq_total_number(gzFile fin)
{
	size_t seq_total_number = 0;
	while (!gzeof(fin)) {
		char buf[256];
		if (!gzgets(fin, buf, sizeof(buf))) break;
		if (buf[0] == '>') {
			++seq_total_number;
		}
	}
	return seq_total_number;
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

	if (output_cmap) {
		size_t seq_total_number = get_seq_total_number(fin);
		gzseek(fin, 0, SEEK_SET);
		write_cmap_header(fout, seq_total_number);
	}

	while (!gzeof(fin)) {
		char buf[256];
		if (!gzgets(fin, buf, sizeof(buf))) break;
		if (buf[0] == '>') {
			print_results(fout, chrom, offset);
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
			offset = 0;
		} else {
			process_line(fout, buf, chrom);
		}
	}
	print_results(fout, chrom, offset);
	gzclose(fout);
	gzclose(fin);
	free_results();
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
