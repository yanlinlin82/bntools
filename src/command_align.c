#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include "nick_map.h"

#define DEF_OUTPUT "stdout"

static int verbose = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools align [options] <input_a> [<input_b>]\n"
			"\n"
			"Options:\n"
			"   <input_x>   input map file(s), in tsv/cmap/bnx format\n"
			"   -v          show verbose message\n"
			"\n");
}

static int align(struct nick_map *map1, struct nick_map *map2)
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
	if (optind >= argc || optind + 2 < argc) {
		print_usage();
		return 1;
	}
	return 0;
}

int align_main(int argc, char * const argv[])
{
	struct nick_map map;
	int ret;

	if (check_options(argc, argv)) {
		return 1;
	}

	nick_map_init(&map);
	if (nick_map_load(&map, argv[optind])) {
		nick_map_free(&map);
		return 1;
	}

	if (optind + 1 >= argc) {
		ret = align(&map, &map);
	} else {
		struct nick_map map2;
		nick_map_init(&map2);
		if (nick_map_load(&map2, argv[optind + 1])) {
			nick_map_free(&map2);
			nick_map_free(&map);
			return 1;
		}
		ret = align(&map, &map2);
		nick_map_free(&map2);
	}
	nick_map_free(&map);
	return ret;
}
