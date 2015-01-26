#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include "nick_map.h"

#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "tsv"

static int verbose = 0;

static char output_file[PATH_MAX] = DEF_OUTPUT;
static int format = FORMAT_TSV;
static char name_list_file[PATH_MAX] = "";

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools view [options] <input> [...]\n"
			"\n"
			"Options:\n"
			"   <input> [...]   input map file(s), in tsv/cmap/bnx/txt format\n"
			"   -o FILE         output file ["DEF_OUTPUT"]\n"
			"   -f STR          output format, tsv/cmap/bnx/txt ["DEF_FORMAT"]\n"
			"   -s FILE         select fragment by name, listed in lines\n"
			"   -v              show verbose message\n"
			"\n");
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "o:f:s:v")) != -1) {
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
		case 's':
			snprintf(name_list_file, sizeof(name_list_file), "%s", optarg);
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
	struct name_list name_list = { };
	int i;

	if (check_options(argc, argv)) {
		return 1;
	}

	if (name_list_file[0]) {
		if (load_name_list(&name_list, name_list_file)) {
			fprintf(stderr, "Error: Can not load name list from '%s'!\n", name_list_file);
			return 1;
		}
	}

	nick_map_init(&map);

	for (i = optind; i < argc; ++i) {
		if (nick_map_load_ex(&map, argv[i],
				(name_list.names.size ? &name_list : NULL))) {
			nick_map_free(&map);
			return 1;
		}
	}

	nick_map_save(&map, output_file, format);
	nick_map_free(&map);
	return 0;
}
