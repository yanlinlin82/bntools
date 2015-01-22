#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "nick_map.h"

#define DEF_OUTPUT "stdout"

static int verbose = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools map [options] <ref> <query>\n"
			"\n"
			"Options:\n"
			"   <ref>     reference genome, in tsv/cmap format\n"
			"   <query>   query molecules/contigs, in tsv/cmap/bnx format\n"
			"   -v        show verbose message\n"
			"\n");
}

static int map(struct nick_map *map1, struct nick_map *map2)
{
	return 0;
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "v")) != -1) {
		switch (c) {
		case 'v':
			++verbose;
			break;
		default:
			return 1;
		}
	}
	if (optind + 2 != argc) {
		print_usage();
		return 1;
	}
	return 0;
}

int map_main(int argc, char * const argv[])
{
	struct nick_map map1, map2;
	int ret;

	if (check_options(argc, argv)) {
		return 1;
	}

	nick_map_init(&map1);
	nick_map_init(&map2);
	if (nick_map_load(&map1, argv[optind])) {
		nick_map_free(&map1);
		return 1;
	}
	if (nick_map_load(&map2, argv[optind + 1])) {
		nick_map_free(&map2);
		nick_map_free(&map1);
		return 1;
	}

	ret = map(&map1, &map2);

	nick_map_free(&map2);
	nick_map_free(&map1);
	return ret;
}
