#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>
#include "nick_map.h"
#include "bn_file.h"

#define DEF_OUTPUT "stdout"
#define DEF_FORMAT "tsv"

struct range {
	char *name;
	int start;
	int end;
};

static int verbose = 0;
static int help = 0;

static char output_file[PATH_MAX] = DEF_OUTPUT;
static int out_format = FORMAT_TSV;
static int counting = 0;
static array(struct range) ranges = { };
static int reverse = 0;

static size_t fragment_count = 0;
static size_t nick_count = 0;
static long long total_size = 0;
static int save_into_map = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools view [options] <input> [...]\n"
			"\n"
			"Options:\n"
			"   <input> [...]  input map file(s), in tsv/cmap/bnx/txt format\n"
			"   -o FILE        output file ["DEF_OUTPUT"]\n"
			"   -f STR         output format, tsv/cmap/bnx/txt ["DEF_FORMAT"]\n"
			"   -r STR         select range(s), specified as string\n"
			"   -R FILE        select range(s), specified as lines in file\n"
			"   -t             transform to reverse order\n"
			"   -c             count fragments, nicks and total size\n"
			"   -v             show verbose message\n"
			"   -h             show this help, '-hh' for more detail help\n"
			"\n");
	if (help > 1) {
		fprintf(stderr, "Note:\n"
				"   Range string is formatted as: <name>:<start>-<end>, where <name>\n"
				"is a string without ':', <start> and <end> are numbers in bp.\n"
				"\n");
	}
}

int ranges_overlap(int start1, int end1, int start2, int end2)
{
	assert(start1 >= 0);
	assert(end1 >= 0);
	assert(start1 <= end1 || end1 == 0);
	assert(start2 >= 0);
	assert(end2 >= 0);
	assert(start2 <= end2 || end2 == 0);

	if (start1 == 0 && end1 == 0) { /* range1 is whole fragment */
		return 1;
	} else if (start2 == 0 && end2 == 0) { /* range2 is whole fragment */
		return 1;
	} else if (start1 == 0) { /* range1 spans to beginning */
		assert(end1 != 0);
		return (end1 > start2);
	} else if (start2 == 0) { /* range2 spans to beginning */
		assert(end2 != 0);
		return (end2 > start1);
	} else if (end1 == 0) { /* range1 spans to end */
		return (start1 < end2 || end2 == 0);
	} else if (end2 == 0) { /* range2 spans to end */
		assert(end1 != 0);
		return (start2 < end1);
	} else {
		assert(start1 != 0 && end1 != 0 && start2 != 0 && end2 != 0);
		return (start1 < end2 && start2 < end1);
	}
}

int append_range(const char *s)
{
	char *name = strdup(s);
	char *p, *q;
	int start = 0, end = 0;
	size_t i;

	p = strchr(name, ':');
	if (p) {
		*p++ = '\0';
		q = strchr(p, '-');
		if (q) {
			*q++ = '\0';
		}
		start = atoi(p);
		end = atoi(q);
	}

	for (i = 0; i < ranges.size; ++i) {
		if (strcmp(ranges.data[i].name, name) == 0) {
			if (ranges_overlap(ranges.data[i].start, ranges.data[i].end, start, end)) {
				fprintf(stderr, "Error: ranges overlap between %s:%d-%d and %s:%d-%d\n",
						name, ranges.data[i].start, ranges.data[i].end,
						name, start, end);
				return 1;
			}
		}
	}

	if (array_reserve(ranges, ranges.size + 1)) {
		return -ENOMEM;
	}
	ranges.data[ranges.size].name = name;
	ranges.data[ranges.size].start = start;
	ranges.data[ranges.size].end = end;
	++ranges.size;
	return 0;
}

int append_ranges_in_file(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!feof(fp)) {
		fprintf(stderr, "Error: Can not open range file '%s'\n", filename);
		return 1;
	}
	while (!feof(fp)) {
		char buf[256];
		if (!fgets(buf, sizeof(buf), fp)) break;
		if (append_range(buf)) {
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);
	return 0;
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "o:f:r:R:tcvh")) != -1) {
		switch (c) {
		case 'o':
			snprintf(output_file, sizeof(output_file), "%s", optarg);
			break;
		case 'f':
			out_format = parse_format_text(optarg);
			if (out_format == FORMAT_UNKNOWN) {
				fprintf(stderr, "Error: Unknown output format '%s'!\n", optarg);
				return 1;
			}
			break;
		case 'r':
			if (append_range(optarg)) {
				return 1;
			}
			break;
		case 'R':
			if (append_ranges_in_file(optarg)) {
				return 1;
			}
			break;
		case 't':
			reverse = 1;
			break;
		case 'c':
			counting = 1;
			break;
		case 'v':
			++verbose;
			break;
		case 'h':
			++help;
			break;
		default:
			return 1;
		}
	}
	if (help || optind >= argc) {
		print_usage();
		return 1;
	}
	return 0;
}

static void free_ranges(void)
{
	size_t i;
	for (i = 0; i < ranges.size; ++i) {
		free(ranges.data[i].name);
	}
	array_free(ranges);
}

static void extract_fragment(const struct fragment *f, int start, int end, struct fragment *sub)
{
	size_t i;

	assert(f != NULL);
	assert(start >= 0);
	assert(end >= 0);
	assert(start <= end || end == 0);
	assert(sub != NULL);

	if (start == 0 && end == 0) {
		snprintf(sub->name, sizeof(sub->name), "%s", f->name);
	} else if (start == 0) {
		snprintf(sub->name, sizeof(sub->name), "%s:-%d", f->name, end);
	} else if (end == 0) {
		snprintf(sub->name, sizeof(sub->name), "%s:%d-", f->name, start);
	} else {
		snprintf(sub->name, sizeof(sub->name), "%s:%d-%d", f->name, start, end);
	}
	sub->size = (end == 0 ? f->size : end) - start;
	for (i = 0; i < f->nicks.size; ++i) {
		if (f->nicks.data[i].pos < start) continue;
		if (end != 0 && f->nicks.data[i].pos > end) break;
		nick_map_add_site(sub, f->nicks.data[i].pos - start, f->nicks.data[i].flag);
	}
}

int process_fragment(struct nick_map *map, struct fragment *f, gzFile file)
{
	if (counting) {
		++fragment_count;
		nick_count += f->nicks.size;
		total_size += f->size;
	} else {
		if (reverse) {
			size_t i;
			int pos;
			for (i = 0; i < f->nicks.size; ++i) {
				f->nicks.data[i].pos = f->size - f->nicks.data[i].pos;
			}
			for (i = 0; i < f->nicks.size / 2; ++i) {
				pos = f->nicks.data[i].pos;
				f->nicks.data[i].pos = f->nicks.data[f->nicks.size - 1 - i].pos;
				f->nicks.data[f->nicks.size - 1 - i].pos = pos;
			}
		}
		if (save_into_map) {
			if (array_reserve(map->fragments, map->fragments.size + 1)) {
				return 1;
			}
			memcpy(&map->fragments.data[map->fragments.size++], f, sizeof(struct fragment));
			f->nicks.data = NULL;
			f->nicks.size = 0;
			f->nicks.capacity = 0;
		} else {
			save_fragment(file, f, out_format);
		}
	}
	return 0;
}

int view_main(int argc, char * const argv[])
{
	struct file *fp;
	struct nick_map map = { };
	struct fragment fragment = { };
	struct fragment sub = { };
	gzFile file;
	int i, j, ret = 0;
	int format;

	if (check_options(argc, argv)) {
		ret = 1;
		goto out;
	}

	if (counting) {
		save_into_map = 0;
		file = NULL;
	} else if (out_format == FORMAT_TXT || out_format == FORMAT_TSV) {
		save_into_map = 0;
		file = open_gzfile_write(output_file);
		if (!file) {
			ret = 1;
			goto out;
		}
		save_header(file, &map, out_format);
	} else {
		save_into_map = 1;
		file = NULL;
	}

	for (i = optind; i < argc; ++i) {
		fp = file_open(argv[i]);
		if (!fp) {
			ret = 1;
			goto out;
		}
		if (bn_read_header(fp, &format, &map) != 0) {
			file_close(fp);
			ret = 1;
			goto out;
		}
		while (bn_read(fp, format, &fragment) == 0) {
			if (ranges.size == 0) {
				if (process_fragment(&map, &fragment, file)) {
					ret = 1;
					goto out;
				}
			} else {
				for (j = 0; j < ranges.size; ++j) {
					if (strcmp(ranges.data[j].name, fragment.name) != 0) continue;
					extract_fragment(&fragment, ranges.data[j].start,
							ranges.data[j].end, &sub);
					if (process_fragment(&map, &sub, file)) {
						ret = 1;
						goto out;
					}
				}
			}
		}
		file_close(fp);
	}

	if (counting) {
		fprintf(stdout, "%zd\t%zd\t%lld\n", fragment_count, nick_count, total_size);
	} else if (save_into_map) {
		nick_map_save(&map, output_file, out_format);
	} else {
		gzclose(file);
	}
out:
	nick_map_free(&map);
	free_ranges();
	return ret;
}
