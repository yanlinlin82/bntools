#include <stdio.h>
#include <unistd.h>
#include "nick_map.h"
#include "ref_map.h"
#include "bn_file.h"

static int verbose = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools motif [options] <ref>\n"
			"\n"
			"Options:\n"
			"   <ref>   reference genome, in tsv/cmap format\n"
			"   -v      show verbose message\n"
			"   -h      show this help\n"
			"\n");
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "vh")) != -1) {
		switch (c) {
		case 'v':
			++verbose;
			break;
		case 'h':
			print_usage();
		default:
			return 1;
		}
	}
	if (optind + 1 != argc) {
		print_usage();
		return 1;
	}
	return 0;
}

int motif_main(int argc, char * const argv[])
{
	struct ref_map ref;

	if (check_options(argc, argv)) {
		return 1;
	}

	ref_map_init(&ref);
	if (nick_map_load(&ref.map, argv[optind])) {
		return 1;
	}
	ref_map_build_index(&ref);
	ref_map_dump(&ref);

	ref_map_free(&ref);
	return 0;
}
