#include <stdio.h>
#include <unistd.h>
#include "nick_map.h"
#include "ref_map.h"

#define DEF_TOLERANCE 0.1
#define DEF_MIN_MATCH 4

static int verbose = 0;
static double tolerance = DEF_TOLERANCE;
static int min_match = DEF_MIN_MATCH;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools map [options] <ref> <query>\n"
			"\n"
			"Options:\n"
			"   <ref>        reference genome, in tsv/cmap format\n"
			"   <query>      query molecules/contigs, in tsv/cmap/bnx format\n"
			"   -v           show verbose message\n"
			"   -e <FLOAT>   tolerance to compare fragment size [%f]\n"
			"   -m <INT>     minimal matched fragment [%d]\n"
			"\n", DEF_TOLERANCE, DEF_MIN_MATCH);
}

static void print_header(void)
{
	printf("#name\tchrom\tpos\tstrand\tsize\tlabels\trstart\trend\tqstart\tqend\n");
}

static void print_item(const char *name, const char *chrom, int pos, int direct,
		int size, size_t labels, size_t rstart, size_t rend, size_t qstart, size_t qend)
{
	printf("%s\t%s\t%d\t%s\t%d\t%zd\t%zd\t%zd\t%zd\t%zd\n",
			name, chrom, pos, (direct > 0 ? "+" : "-"),
			size, labels, rstart, rend, qstart, qend);
}

static void map(struct ref_map *ref, struct fragment *qry_item)
{
	size_t i, j, index;
	int *intervals;
	int direct;
	size_t qstart = 1;

	if (qry_item->nicks.size < min_match + 2) {
		if (verbose) {
			fprintf(stderr, "Warning: Skip '%s' for less than %d labels!\n", qry_item->name, min_match);
		}
		return;
	}

	intervals = malloc(sizeof(int) * qry_item->nicks.size);
	for (i = 0; i < qry_item->nicks.size; ++i) {
		intervals[i] = qry_item->nicks.data[i].pos - (i == 0 ? 0 : qry_item->nicks.data[i - 1].pos);
	}

	while (qstart + 1 < qry_item->nicks.size) {
		size_t max_count = 0;
		for (i = 0; i < ref->size; ++i) {
			if (ref->index[i]->size < intervals[qstart] * (1 - tolerance)) continue;
			if (ref->index[i]->size > intervals[qstart] * (1 + tolerance)) break;

			index = ref->index[i] - ref->nodes;

			for (direct = 1; direct >= -1; direct -= 2) {
				for (j = 1; j < qry_item->nicks.size; ++j) {
					if (ref->nodes[index + direct * j].size < intervals[qstart + j] * (1 - tolerance)) break;
					if (ref->nodes[index + direct * j].size > intervals[qstart + j] * (1 + tolerance)) break;
				}
				if (j >= min_match) {
					int size = ref->nodes[index + direct * j].pos - ref->nodes[index].pos;
					if (size < 0) size = -size;

					print_item(qry_item->name,
							ref->map.fragments.data[ref->index[i]->chrom].name,
							ref->index[i]->pos, direct,
							size, j, index, index + direct * (j - 1),
							qstart, qstart + j - 1);
					if (max_count < j) {
						max_count = j;
					}
				}
			}
		}
		qstart += (max_count > 0 ? max_count : 1);
	}
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "ve:m:")) != -1) {
		switch (c) {
		case 'v':
			++verbose;
			break;
		case 'e':
			tolerance = atof(optarg);
			break;
		case 'm':
			min_match = atoi(optarg);
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
	struct ref_map ref;
	struct nick_map qry;
	size_t i;

	if (check_options(argc, argv)) {
		return 1;
	}

	ref_map_init(&ref);
	if (nick_map_load(&ref.map, argv[optind])) {
		return 1;
	}
	ref_map_build_index(&ref);

	nick_map_init(&qry);
	if (nick_map_load(&qry, argv[optind + 1])) {
		ref_map_free(&ref);
		return 1;
	}

	print_header();
	for (i = 0; i < qry.fragments.size; ++i) {
		map(&ref, &qry.fragments.data[i]);
	}

	nick_map_free(&qry);
	ref_map_free(&ref);
	return 0;
}
