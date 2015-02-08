#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>
#include "base_map.h"
#include "ref_map.h"
#include "bn_file.h"

#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "tsv"
#define DEF_ENZ_NAME "BspQI"
#define DEF_REC_SEQ "GCTCTTCN^"

static int verbose = 0;

static char enzyme[MAX_ENZYME_NAME_SIZE] = DEF_ENZ_NAME;
static char rec_seq[MAX_REC_SEQ_SIZE] = DEF_REC_SEQ;
static char output_file[PATH_MAX] = DEF_OUTPUT;
static int format = FORMAT_TSV;
static int chrom_only = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools nick [options] <x.fa> [...]\n"
			"\n"
			"Options:\n"
			"   <x.fa> [...]   input FASTA file(s) to generate restriction map\n"
			"   -o FILE        output file ["DEF_OUTPUT"]\n"
			"   -f STR         output format, tsv/cmap/bnx/txt ["DEF_FORMAT"]\n"
			"   -e STR         restriction enzyme name ["DEF_ENZ_NAME"]\n"
			"   -r STR         recognition sequence ["DEF_REC_SEQ"]\n"
			"   -S             select only chr1-22, chrX and chrY to nick\n"
			"   -v             show verbose messages\n"
			"   -h             show this help\n"
			"\n");
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "o:f:e:r:Svh")) != -1) {
		switch (c) {
		case 'o':
			snprintf(output_file, sizeof(output_file), "%s", optarg);
			break;
		case 'f':
			format = parse_format_text(optarg);
			if (format == FORMAT_UNKNOWN) {
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
		case 'S':
			chrom_only = 1;
			break;
		case 'v':
			++verbose;
			break;
		case 'h':
			print_usage();
		default:
			return 1;
		}
	}
	if (optind >= argc) {
		print_usage();
		return 1;
	}
	return 0;
}

int nick_main(int argc, char * const argv[])
{
	struct rec_site site;
	struct ref_map ref;
	int i, ret = 0;

	if (check_options(argc, argv)) {
		return 1;
	}

	ref_map_init(&ref);

	if ((ret = prepare_rec_site(&site, enzyme, rec_seq)) != 0) {
		goto final;
	}

	for (i = optind; i < argc; ++i) {
		if ((ret = nick_map_load_fasta(&ref, argv[i],
				&site, chrom_only, verbose)) != 0) {
			goto final;
		}
	}
	ret = nick_map_save(&ref.map, output_file, format);
final:
	ref_map_free(&ref);
	return ret;
}
