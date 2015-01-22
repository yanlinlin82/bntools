#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "nick_map.h"

#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "tsv"

static int verbose = 0;

static char output_file[PATH_MAX] = DEF_OUTPUT;
static int output_cmap = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools view [options] <input> [...]\n"
			"\n"
			"Options:\n"
			"   <input> [...]   input map file(s), in tsv/cmap/bnx format\n"
			"   -o FILE         output file ["DEF_OUTPUT"]\n"
			"   -f {tsv|cmap}   output format ["DEF_FORMAT"]\n"
			"   -v              show verbose message\n"
			"\n");
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "o:f:v")) != -1) {
		switch (c) {
		case 'o':
			snprintf(output_file, sizeof(output_file), "%s", optarg);
			break;
		case 'f':
			if (strcmp(optarg, "tsv") == 0) {
				output_cmap = 0;
			} else if (strcmp(optarg, "cmap") == 0) {
				output_cmap = 1;
			} else {
				fprintf(stderr, "Error: Unknown output format '%s'!\n", optarg);
				return 1;
			}
			break;
		case 'v':
			++verbose;
			break;
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

int view_main(int argc, char * const argv[])
{
	struct nick_map map;
	int i;

	if (check_options(argc, argv)) {
		return 1;
	}

	nick_map_init(&map);

	for (i = optind; i < argc; ++i) {
		if (nick_map_load(&map, argv[i])) {
			nick_map_free(&map);
			return 1;
		}
	}

	nick_map_save(&map, output_file, output_cmap);
	nick_map_free(&map);
	return 0;
}
