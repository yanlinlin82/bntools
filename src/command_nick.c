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

static char enzyme[MAX_ENZYME_NAME_SIZE] = DEF_ENZ_NAME;
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
			char c = char_to_base(rec_seq[i]);
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
		if (rec_bases[i] != base_to_comp(rec_bases[rec_seq_size - i - 1])) {
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

static int seq_match(const char *ref, const char *query, size_t len, int strand)
{
	size_t i;
	assert(strand == STRAND_PLUS || strand == STRAND_MINUS);
	for (i = 0; i < len; ++i) {
		char r = ref[i];
		char q = (strand == STRAND_PLUS ? query[i] : base_to_comp(query[len - i - 1]));
		if ((r & q) != r) {
			return 0;
		}
	}
	return 1;
}

static int process_line(struct nick_list *list, const char *line, const char *chrom, int base_count)
{
	const char *p;
	int strand;
	int matched;
	for (p = line; *p; ++p) {
		if (is_whitespace(*p)) {
			continue;
		}
		if (pos >= sizeof(buf)) {
			memcpy(buf, buf + sizeof(buf) - rec_seq_size + 1, rec_seq_size - 1);
			pos = rec_seq_size - 1;
		}
		++base_count;
		buf[pos++] = char_to_base(*p);
		if (pos < rec_seq_size) {
			continue;
		}
		for (strand = STRAND_PLUS, matched = 0; strand <= STRAND_MINUS; ++strand) {
			if (matched || seq_match(buf + pos - rec_seq_size, rec_bases, rec_seq_size, strand)) {
				int site_pos = base_count - (strand == STRAND_MINUS ? nick_offset : (rec_seq_size - nick_offset));
				if (nick_map_add_site(list, site_pos, strand)) {
					return -ENOMEM;
				}
				matched = palindrome;
			}
		}
	}
	return base_count;
}

static int process(gzFile fin, gzFile fout)
{
	struct nick_list *list = NULL;
	char chrom[MAX_CHROM_NAME_SIZE] = "";
	int base_count = 0;
	while (!gzeof(fin)) {
		char buf[256];
		if (!gzgets(fin, buf, sizeof(buf))) break;
		if (buf[0] == '>') {
			if (list) {
				if (!output_cmap) {
					nick_map_write(fout, &map, list);
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
			nick_map_write(fout, &map, list);
		}
	}
	if (output_cmap) {
		nick_map_write_cmap(fout, &map);
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
			snprintf(enzyme, sizeof(enzyme), "%s", optarg);
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

	nick_map_init(&map, enzyme, rec_seq);

	return nick(argv[optind]);
}
