#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <time.h>
#include <zlib.h>
#include "nick_map.h"

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

#define FIRST_FRAGMENT 1
#define LAST_FRAGMENT  2

struct node {
	size_t chrom;  /* item index in ref */
	int pos;
	int size;
	int flag;
};

static size_t ref_node_count = 0;
static struct node * ref_nodes = NULL;
static struct node **ref_index = NULL;

int sort_by_size(const void *a, const void *b)
{
	size_t i, j, k;
	const struct node * const *node_a = a;
	const struct node * const *node_b = b;
	if ((*node_a)->size < (*node_b)->size) return -1;
	if ((*node_a)->size > (*node_b)->size) return 1;
	i = (*node_a) - ref_nodes;
	j = (*node_b) - ref_nodes;
	for (k = 0; ; ++k) {
		if (ref_nodes[i + k].size < ref_nodes[j + k].size) return -1;
		if (ref_nodes[i + k].size > ref_nodes[j + k].size) return 1;
	}
	return 0;
}

static void generate_ref_nodes(const struct nick_map *ref)
{
	size_t i, j, k;

	ref_node_count = 0;
	for (i = 0; i < ref->size; ++i) {
		ref_node_count += ref->data[i].size;
	}

	ref_nodes = malloc(sizeof(struct node) * ref_node_count);
	ref_index = malloc(sizeof(struct node *) * ref_node_count);

	for (i = 0, k = 0; i < ref->size; ++i) {
		const struct nick_list *list = &ref->data[i];
		for (j = 0; j < list->size; ++j) {
			ref_nodes[k].chrom = i;
			ref_nodes[k].pos = list->data[j].pos;
			ref_nodes[k].size = list->data[j].pos - (j == 0 ? 0 : list->data[j - 1].pos);
			ref_nodes[k].flag = (j == 0 ? FIRST_FRAGMENT : 0) | (j + 1 == list->size ? LAST_FRAGMENT : 0);
			ref_index[k] = &ref_nodes[k];
			++k;
		}
	}

	qsort(ref_index, ref_node_count, sizeof(struct node *), sort_by_size);

	if (verbose > 0) {
		printf("Total fragment count = %zd\n", ref_node_count);
		for (i = 0; i < ref_node_count; ++i) {
			printf("  [%2zd] %s, %d, %d, %d\n", i,
					ref->data[ref_index[i]->chrom].fragment_name,
					ref_index[i]->pos, ref_index[i]->size, ref_index[i]->flag);
		}
	}
}

#define TOLERANCE 0.1
#define MIN_MATCH 4

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

static void map(struct nick_map *ref, struct nick_list *qry_item)
{
	size_t i, j, index;
	int *intervals;
	int direct;
	size_t qstart = 1;

	if (qry_item->size < MIN_MATCH) {
		if (verbose) {
			fprintf(stderr, "Warning: Skip '%s' for less than %d labels!\n", qry_item->fragment_name, MIN_MATCH);
		}
		return;
	}

	intervals = malloc(sizeof(int) * qry_item->size);
	for (i = 0; i < qry_item->size; ++i) {
		intervals[i] = qry_item->data[i].pos - (i == 0 ? 0 : qry_item->data[i - 1].pos);
	}

	while (qstart + 1 < qry_item->size) {
		size_t max_count = 0;
		for (i = 0; i < ref_node_count; ++i) {
			if (ref_index[i]->size < intervals[qstart] * (1 - TOLERANCE)) continue;
			if (ref_index[i]->size > intervals[qstart] * (1 + TOLERANCE)) break;

			index = ref_index[i] - ref_nodes;

			for (direct = 1; direct >= -1; direct -= 2) {
				for (j = 1; j < qry_item->size; ++j) {
					if (ref_nodes[index + direct * j].size < intervals[qstart + j] * (1 - TOLERANCE)) break;
					if (ref_nodes[index + direct * j].size > intervals[qstart + j] * (1 + TOLERANCE)) break;
				}
				if (j >= MIN_MATCH) {
					int size = ref_nodes[index + direct * j].pos - ref_nodes[index].pos;
					if (size < 0) size = -size;

					print_item(qry_item->fragment_name,
							ref->data[ref_index[i]->chrom].fragment_name,
							ref_index[i]->pos, direct,
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
	struct nick_map ref, qry;
	size_t i;

	if (check_options(argc, argv)) {
		return 1;
	}

	nick_map_init(&ref);
	if (nick_map_load(&ref, argv[optind])) {
		return 1;
	}
	{
		clock_t t0 = clock();
		generate_ref_nodes(&ref);
		printf("generating ref_nodes used %ld ms.\n",
				(clock() - t0) * 1000 / CLOCKS_PER_SEC);
	}

	nick_map_init(&qry);
	if (nick_map_load(&qry, argv[optind + 1])) {
		nick_map_free(&ref);
		return 1;
	}

	print_header();
	for (i = 0; i < qry.size; ++i) {
		map(&ref, &qry.data[i]);
	}

	nick_map_free(&qry);
	nick_map_free(&ref);
	return 0;
}
