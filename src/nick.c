#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "nick.h"

#define NAME "nick"

static int verbose = 0;

static char enzyme_name[16] = "BspQI";
static char recog_seq[32] = "GCTCTTCN^";
static char output_file[PATH_MAX] = "stdout";
static int output_cmap = 0;

static char target_seq[32] = "";
static int target_length = 0;
static int target_offset = -1;

static char target_seq_revcomp[32] = "";
static int target_offset_revcomp = -1;

static int chrom = 0;
static int offset = 0;

static char buf[64];
static int pos = 0;

static int prepare_recog_seq(void)
{
	int i;
	for (i = 0, target_length = 0; recog_seq[i]; ++i) {
		if (recog_seq[i] == '^') {
			if (target_offset >= 0) {
				fprintf(stderr, NAME": invalid recognition sequence '%s'\n", recog_seq);
				return 1;
			}
			target_offset = i;
		} else {
			if (target_length >= sizeof(target_seq)) {
				fprintf(stderr, NAME": recognition sequence is too long '%s'\n", recog_seq); 
				return 1;
			}
			target_seq[target_length++] = recog_seq[i];
		}
	}
	target_seq[target_length] = '\0';

	for (i = 0; i < target_length; ++i) {
		switch (target_seq[i]) {
		case 'A': case 'a': target_seq_revcomp[target_length - i - 1] = 'T'; break;
		case 'C': case 'c': target_seq_revcomp[target_length - i - 1] = 'G'; break;
		case 'G': case 'g': target_seq_revcomp[target_length - i - 1] = 'C'; break;
		case 'T': case 't': target_seq_revcomp[target_length - i - 1] = 'A'; break;
		default: target_seq_revcomp[target_length - i - 1] = 'N'; break;
		}
	}
	if (memcmp(target_seq, target_seq_revcomp, target_length) != 0) {
		target_offset_revcomp = target_length - target_offset;
	}
	return 0;
}

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools "NAME" [options] <fa.gz>\n"
			"\n"
			"Options:\n"
			"   -e STR   restriction enzyme name [BspQI]\n"
			"   -r STR   recognition sequence [GCTCTTCN^]\n"
			"   -o FILE  output file [stdout]\n"
			"   -C       output in .cmap format rather than .bed format\n"
			"   -v       show verbose message\n"
			"\n");
}

static inline int is_whitespace(int c)
{
	return (c == ' ' || c == '\t' || c == '\r' || c == '\n');
}

static int load_name(const char *buf, char **end, char *name, size_t name_len)
{
	const char *p = buf;
	int chrom = 0;

	if (name && name_len > 0) {
		size_t i;
		for (i = 0; i < name_len; ++i) {
			if (buf[i] == ' ' || buf[i] == '\t' || buf[i] == '\n' || buf[i] == '\r') {
				break;
			}
			name[i] = buf[i];
		}
		if (i >= name_len) {
			i = name_len - 1;
		}
		name[i] = '\0';
	}

	if (memcmp(p, "chr", 3) == 0) {
		p += 3;
	}
	if (*p >= '1' && *p <= '9') {
		chrom = (*p++ - '0');
		while (*p >= '0' && *p <= '9') {
			chrom = (chrom * 10) + (*p++ - '0');
		}
		if (chrom > 22 || !is_whitespace(*p)) {
			return 0;
		}
	} else if ((*p == 'X' || *p == 'x') && is_whitespace(p[1])) {
		++p;
		chrom = 23;
	} else if ((*p == 'Y' || *p == 'y') && is_whitespace(p[1])) {
		++p;
		chrom = 24;
	} else {
		return 0;
	}
	if (end) {
		*end = (char *)p;
	}
	return chrom;
}

struct integer_array {
	size_t capacity;
	size_t size;
	int *data;
};

static struct integer_array sites = { 0 };

static void free_sites(void)
{
	free(sites.data);
	sites.data = NULL;
	sites.size = 0;
	sites.capacity = 0;
}

static size_t new_capacity(size_t capacity, size_t size)
{
	while (capacity < size) {
		if (capacity < 16) {
			capacity += 16;
		} else if (capacity > 64) {
			capacity += 64;
		} else {
			capacity += capacity;
		}
	}
	return capacity;
}

static int append_sites(int site)
{
	if (sites.size == sites.capacity) {
		size_t capacity = new_capacity(sites.capacity, sites.size + 1);
		int *p = (int *)malloc(sizeof(int) * capacity);
		if (!p) {
			return -ENOMEM;
		}
		if (sites.size > 0) {
			memcpy(p, sites.data, sizeof(int) * sites.size);
		}
		free(sites.data);
		sites.data = p;
		sites.capacity = capacity;
	}
	sites.data[sites.size++] = site;
	return 0;
}

static int seq_match(const char *ref, const char *query, size_t len)
{
	size_t i;
	for (i = 0; i < len; ++i) {
		if (query[i] == 'N' || query[i] == 'n') continue;
		if ((ref[i] == 'A' || ref[i] == 'a') && (query[i] == 'A' || query[i] == 'a')) continue;
		if ((ref[i] == 'C' || ref[i] == 'c') && (query[i] == 'C' || query[i] == 'c')) continue;
		if ((ref[i] == 'G' || ref[i] == 'g') && (query[i] == 'G' || query[i] == 'g')) continue;
		if ((ref[i] == 'T' || ref[i] == 't') && (query[i] == 'T' || query[i] == 't')) continue;
		return 0;
	}
	return 1;
}

static void process_line(gzFile fout, const char *line, const char *chrom)
{
	const char *p;
	for (p = line; *p; ++p) {
		if (is_whitespace(*p)) {
			continue;
		}
		if (pos >= sizeof(buf)) {
			memcpy(buf, buf + sizeof(buf) - target_length + 1, target_length - 1);
			pos = target_length - 1;
		}
		++offset;
		buf[pos++] = *p;
		if (pos >= target_length) {
			if (seq_match(buf + pos - target_length, target_seq, target_length)) {
				int site = offset - target_length + target_offset;
				if (output_cmap) {
					append_sites(site);
				} else {
					gzprintf(fout, "%s\t%d\t%d\t%s\t0\t+\n", chrom, site - 1, site, enzyme_name);
				}
			} else {
				if (target_offset_revcomp >= 0) {
					if (seq_match(buf + pos - target_length, target_seq_revcomp, target_length)) {
						int site = offset - target_length + target_offset_revcomp;
						if (output_cmap) {
							append_sites(offset - target_length + target_offset_revcomp);
						} else {
							gzprintf(fout, "%s\t%d\t%d\t%s\t0\t-\n", chrom, site - 1, site, enzyme_name);
						}
					}
				}
			}
		}
	}
}

static void write_cmap_header(gzFile fout)
{
	gzprintf(fout, "# CMAP File Version:  0.1\n");
	gzprintf(fout, "# Label Channels:  1\n");
	gzprintf(fout, "# Nickase Recognition Site 1:  %s/%s\n", enzyme_name, recog_seq);
	gzprintf(fout, "# Number of Consensus Nanomaps:    24\n");
	gzprintf(fout, "#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence\n");
	gzprintf(fout, "#f int	float	int	int	int	float	float	int	int\n");
}

static void output_cmap_line(gzFile fout, int cmap_id, int contig_length, size_t num_sites, size_t site_id,
		int label_channel, int position, int stddev, int coverage, int occurance)
{
	gzprintf(fout, "%d\t%d\t%zd\t%zd\t%d\t%d\t%d\t%d\t%d\n",
			cmap_id, contig_length, num_sites, site_id,
			label_channel, position, stddev, coverage, occurance);
}

static void output_chrom_sites(gzFile fout)
{
	size_t i;
	if (verbose > 0) {
		fprintf(stderr, "%d bp\n", offset);
	}
	for (i = 0; i < sites.size; ++i) {
		output_cmap_line(fout, chrom, offset, sites.size, i + 1, 1, sites.data[i], 0, 0, 0);
	}
	output_cmap_line(fout, chrom, offset, sites.size, sites.size + 1, 0, offset, 0, 1, 1);
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
		fprintf(stderr, NAME": Can not open FASTA file '%s'\n", in);
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
			fprintf(stderr, NAME": Output file '%s' has already existed!\n", output_file);
		} else {
			fprintf(stderr, NAME": Can not open output file '%s'\n", output_file);
		}
		gzclose(fin);
		return 1;
	}

	c = gzgetc(fin);
	if (c != '>') {
		fprintf(stderr, NAME": File '%s' is not in FASTA format\n", in);
		gzclose(fout);
		gzclose(fin);
		return 1;
	}
	gzungetc(c, fin);

	if (output_cmap) {
		write_cmap_header(fout);
	}

	char name[64] = "";
	while (!gzeof(fin)) {
		char buf[256];
		if (!gzgets(fin, buf, sizeof(buf))) break;
		if (buf[0] == '>') {
			char *p = NULL;
			if (chrom > 0 && output_cmap) {
				output_chrom_sites(fout);
			}
			chrom = load_name(buf + 1, &p, name, sizeof(name));
			if (chrom > 0) {
				*p = '\0';
				if (verbose > 0) {
					fprintf(stderr, "Loading sequence '%s' ... ", buf + 1);
				}
				offset = 0;
				sites.size = 0;
			}
		} else {
			process_line(fout, buf, name);
		}
	}
	if (chrom > 0 && output_cmap) {
		output_chrom_sites(fout);
	}
	gzclose(fout);
	gzclose(fin);
	free_sites();
	return 0;
}

int nick_main(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "e:r:v:o:C")) != -1) {
		switch (c) {
		case 'e':
			snprintf(enzyme_name, sizeof(enzyme_name), "%s", optarg);
			break;
		case 'r':
			snprintf(recog_seq, sizeof(recog_seq), "%s", optarg);
			break;
		case 'o':
			snprintf(output_file, sizeof(output_file), "%s", optarg);
			break;
		case 'v':
			++verbose;
			break;
		case 'C':
			output_cmap = 1;
			break;
		default:
			return 1;
		}
	}
	if (argc != optind + 1) {
		print_usage();
		return 1;
	}
	if (prepare_recog_seq() != 0) {
		return 1;
	}
	return nick(argv[optind]);
}
