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
			"Usage: bntools map [options] <ref.fa> <mol.bnx>\n"
			"\n"
			"Options:\n"
			"   <ref.fa>   reference genome, in fasta/tsv/cmap format\n"
			"   <mol.bnx>  bionano molecules, in bnx/tsv/cmap format\n"
			"   -o FILE    output file [stdout]\n"
			"   -v         show verbose message\n"
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

static int map(struct nick_map *map1, struct nick_map *map2)
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

	gzclose(out);
	return 0;
}

static int process(const char *ref, const char *mol)
{
	struct nick_map map1, map2;
	int ret;

	nick_map_init(&map1);
	nick_map_init(&map2);

	if (load_map(&map1, ref)) {
		nick_map_free(&map1);
		return 1;
	}
	if (load_map(&map2, mol)) {
		nick_map_free(&map2);
		nick_map_free(&map1);
		return 1;
	}

	ret = map(&map1, &map2);

	nick_map_free(&map2);
	nick_map_free(&map1);

	return ret;
}

int map_main(int argc, char * const argv[])
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
	if (optind + 2 != argc) {
		print_usage();
		return 1;
	}
	return process(argv[optind], argv[optind + 1]);
}
