#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>
#include "nick_map.h"
#include "bn_file.h"

#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "tsv"

static int verbose = 0;

static char output_file[PATH_MAX] = DEF_OUTPUT;
static int format = FORMAT_TSV;
static char name_list_file[PATH_MAX] = "";
static int counting = 0;

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
			"   -c              count fragments and nicks\n"
			"   -v              show verbose message\n"
			"   -h              show this help\n"
			"\n");
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "o:f:s:cvh")) != -1) {
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
		case 'c':
			counting = 1;
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

int view_main(int argc, char * const argv[])
{
	struct bn_file *fp;
	struct nick_map map;
	struct fragment fragment = { };
	struct name_list name_list = { };
	size_t fragment_count = 0;
	size_t nick_count = 0;
	gzFile file;
	int i, ret = 0, save_into_map;

	if (check_options(argc, argv)) {
		return 1;
	}

	nick_map_init(&map);

	if (counting) {
		save_into_map = 0;
		file = NULL;
	} else if (format == FORMAT_TXT || format == FORMAT_TSV) {
		save_into_map = 0;
		file = open_gzfile_write(output_file);
		if (!file) {
			return 1;
		}
		save_header(file, &map, format);
	} else {
		save_into_map = 1;
		file = NULL;
	}

	if (name_list_file[0]) {
		if (load_name_list(&name_list, name_list_file)) {
			fprintf(stderr, "Error: Can not load name list from '%s'!\n", name_list_file);
			return 1;
		}
	}

	for (i = optind; i < argc; ++i) {
		fp = bn_open(argv[i]);
		if (!fp) {
			ret = 1;
			goto out;
		}
		while (bn_read(fp, &fragment) == 0) {
			if (name_list_file[0]) {
				if (!name_list_has(&name_list, fragment.name)) {
					continue;
				}
			}
			if (counting) {
				++fragment_count;
				nick_count += fragment.nicks.size;
			} else if (save_into_map) {
				if (array_reserve(map.fragments, map.fragments.size + 1)) {
					nick_map_free(&map);
					return -ENOMEM;
				}
				memcpy(&map.fragments.data[map.fragments.size++],
						&fragment, sizeof(struct fragment));
			} else {
				save_fragment(file, &fragment, format);
			}
		}
		bn_close(fp);
	}

	if (counting) {
		fprintf(stdout, "%zd\t%zd\n", fragment_count, nick_count);
	} else if (save_into_map) {
		nick_map_save(&map, output_file, format);
	} else {
		gzclose(file);
	}
out:
	nick_map_free(&map);
	return ret;
}
