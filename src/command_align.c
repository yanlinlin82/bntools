#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "nick_map.h"

static int verbose = 0;

static char output_file[PATH_MAX] = "stdout";

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools align [options] <1.map> [<2.map>]\n"
			"\n"
			"Options:\n"
			"   <x.map>  input map file, in tsv/cmap format\n"
			"   -o FILE  output file [stdout]\n"
			"   -v       show verbose message\n"
			"\n");
}

static int load_map(struct nick_map *map, const char *filename)
{
	gzFile file;
	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		file = gzdopen(0, "r"); /* stdin */
	} else {
		file = gzopen(filename, "r");
	}
	if (!file) {
		fprintf(stderr, "Error: Can not open FASTA file '%s'\n", filename);
		return 1;
	}

	if (nick_map_load(map, file)) {
		gzclose(file);
		return 1;
	}

	gzclose(file);
	return 0;
}

static int align(struct nick_map *map1, struct nick_map *map2)
{
	gzFile out;
	if (strcmp(output_file, "-") == 0 || strcmp(output_file, "stdout") == 0) {
		out = gzdopen(1, "wT"); /* stdout, without compression */
	} else {
		size_t len = strlen(output_file);
		if (len > 3 && strcmp(output_file + len - 3, ".gz") == 0) {
			out = gzopen(output_file, "wx"); /* 'x' is for checking existance */
		} else {
			out = gzopen(output_file, "wxT"); /* without compression */
		}
	}
	if (!out) {
		if (errno == EEXIST) {
			fprintf(stderr, "Error: Output file '%s' has already existed!\n", output_file);
		} else {
			fprintf(stderr, "Error: Can not open output file '%s'\n", output_file);
		}
		return 1;
	}

	nick_map_write_cmap(out, map1);

	gzclose(out);
	return 0;
}

int align_main(int argc, char * const argv[])
{
	int c, ret;
	while ((c = getopt(argc, argv, "o:v")) != -1) {
		switch (c) {
		case 'o':
			snprintf(output_file, sizeof(output_file), "%s", optarg);
			break;
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

	{
		struct nick_map map;
		nick_map_init(&map);
		if (load_map(&map, argv[optind])) {
			nick_map_free(&map);
			return 1;
		}
		if (optind + 1 >= argc) {
			ret = align(&map, &map);
		} else {
			struct nick_map map2;
			nick_map_init(&map2);
			if (load_map(&map2, argv[optind + 1])) {
				nick_map_free(&map2);
				nick_map_free(&map);
				return 1;
			}
			ret = align(&map, &map2);
			nick_map_free(&map2);
		}
		nick_map_free(&map);
	}
	return ret;
}
