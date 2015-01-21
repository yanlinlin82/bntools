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
			"Usage: bntools align [options] <1.map> <2.map>\n"
			"\n"
			"Options:\n"
			"   <x.map>  input map file, in tsv/cmap format\n"
			"   -o FILE  output file [stdout]\n"
			"   -v       show verbose message\n"
			"\n");
}

static int align(const char *map1, const char *map2)
{
	gzFile inMap1;
	gzFile inMap2;
	gzFile out;
	struct nick_map map;

	if (strcmp(map1, "-") == 0 || strcmp(map1, "stdin") == 0) {
		inMap1 = gzdopen(0, "r"); /* stdin */
	} else {
		inMap1 = gzopen(map1, "r");
	}
	if (!inMap1) {
		fprintf(stderr, "Error: Can not open FASTA file '%s'\n", map1);
		return 1;
	}

	if (strcmp(map2, "-") == 0 || strcmp(map2, "stdin") == 0) {
		inMap2 = gzdopen(0, "r"); /* stdin */
	} else {
		inMap2 = gzopen(map2, "r");
	}
	if (!inMap2) {
		fprintf(stderr, "Error: Can not open FASTA file '%s'\n", map2);
		gzclose(inMap1);
		return 1;
	}

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
		gzclose(inMap2);
		gzclose(inMap1);
		return 1;
	}

	nick_map_init(&map);

	if (nick_map_load(&map, inMap1)) {
		gzclose(out);
		gzclose(inMap2);
		gzclose(inMap1);
		return 1;
	}

	nick_map_write_cmap(out, &map);

	gzclose(out);
	gzclose(inMap2);
	gzclose(inMap1);
	return 0;
}

int align_main(int argc, char * const argv[])
{
	int c;
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
	if (argc != optind + 2) {
		print_usage();
		return 1;
	}
	return align(argv[optind], argv[optind + 1]);
}
