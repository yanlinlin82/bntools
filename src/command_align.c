#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include "nick_map.h"

#define DEF_OUTPUT "stdout"

static int verbose = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools align [options] <map_a> [<map_b>]\n"
			"\n"
			"Options:\n"
			"   <map_a/b>   input map file(s), in tsv/cmap/bnx format\n"
			"   -v          show verbose message\n"
			"\n");
}

#define TOLERANCE .1
#define MIN_SCORE 2

struct node {
	int score;
	int delta_i;
	int delta_j;
};

static inline int min(int a, int b) { return (a < b ? a : b); }
static inline int max(int a, int b) { return (a > b ? a : b); }

static int align(const char *name_a, const int *a, size_t size_a,
		const char *name_b, const int *b, size_t size_b)
{
	struct node *m;
	size_t size, i, j, i_start, j_start;
	size_t *result_a, *result_b, result_a_count, result_b_count;

	size = sizeof(struct node) * size_a * size_b;
	m = malloc(size);
	if (!m) {
		fprintf(stderr, "Error: Failed to allocate memory for DP matrix!\n");
		return -ENOMEM;
	}
	memset(m, 0, size);

	if (verbose) {
		fprintf(stdout, "align between '%s' and '%s'\n", name_a, name_b);
	}

	for (j = 1; j < size_b; ++j) {
		for (i = 1; i < size_a; ++i) {
			size_t index = j * size_a + i;

			for (j_start = j; j_start > 0; --j_start) {
				int frag_b = b[j] - b[j_start - 1];

				for (i_start = i; i_start > 0; --i_start) {
					int frag_a = a[i] - a[i_start - 1];

					if (j_start + 2 < j || i_start + 2 < i) continue;
					//if (j_start + 1 < j && i_start + 1 < i) continue;
					if (frag_b > frag_a * (1 + TOLERANCE)) break;

					if (abs(frag_a - frag_b) < min(frag_a, frag_b) * TOLERANCE) {
						int score = m[(j_start - 1) * size_a + (i_start - 1)].score + 1;
						if (m[index].score < score) {
							m[index].score = score;
							m[index].delta_i = i - i_start + 1;
							m[index].delta_j = j - j_start + 1;
						}
					}
				}
			}
		}
	}

	if (verbose > 0) {
		printf("a (%zd):\n", size_a);
		for (i = 1; i < size_a; ++i) {
			printf(" [%2zd] %-5d", i, a[i] - a[i - 1]);
			if (i % 5 == 0) printf("\n");
		}
		if (i % 5 != 1) printf("\n");

		printf("b (%zd):\n", size_b);
		for (i = 1; i < size_b; ++i) {
			printf(" [%2zd] %-5d", i, b[i] - b[i - 1]);
			if (i % 5 == 0) printf("\n");
		}
		if (i % 5 != 1) printf("\n");

		printf("matrix:\n");
		for (j = 0; j < size_b; ++j) {
			for (i = 0; i < size_a; ++i) {
				const struct node *p = &m[j * size_a + i];
				printf("%2d/%d,%d", p->score, p->delta_i, p->delta_j);
			}
			printf("\n");
		}
	}

	result_a = malloc(sizeof(size_t) * size_a);
	result_b = malloc(sizeof(size_t) * size_b);
	for (;;) {
		int max_score = 0;
		size_t max_score_i = 0;
		size_t max_score_j = 0;

		for (j = 0; j < size_b; ++j) {
			for (i = 0; i < size_a; ++i) {
				if (max_score < m[j * size_a + i].score) {
					max_score = m[j * size_a + i].score;
					max_score_i = i;
					max_score_j = j;
				}
			}
		}

		if (max_score < MIN_SCORE) break;

		if (verbose > 0) {
			printf("---- max: %d (%zd, %zd)\n", max_score, max_score_i, max_score_j);
		}
		i = max_score_i;
		j = max_score_j;
		result_a_count = 0;
		result_b_count = 0;
		for (;;) {
			size_t k;
			struct node *p = &m[j * size_a + i];

			result_a[result_a_count++] = i;
			result_b[result_b_count++] = j;

			if (p->delta_i == 0 && p->delta_j == 0) break;

			if (verbose > 0) {
				printf("a: %d { ", a[i] - a[i - p->delta_i]);
				for (k = 0; k < p->delta_i; ++k) {
					size_t index = i - p->delta_i + k;
					if (k > 0) printf(", ");
					printf("[%zd] %d", index, a[index + 1] - a[index]);
				}
				printf(" }\t");

				printf("b: %d { ", b[i] - b[i - p->delta_j]);
				for (k = 0; k < p->delta_j; ++k) {
					size_t index = j - p->delta_j + k;
					if (k > 0) printf(", ");
					printf("[%zd] %d", index, b[index + 1] - b[index]);
				}
				printf(" }\n");
			}

			i -= p->delta_i;
			j -= p->delta_j;
			p->score = 0;
		}

		{
			static int count = 0;
			fprintf(stdout, ">0\t%d\t%s\t%s\n", ++count, name_a, name_b);
		}
		for (i = 0; i < result_a_count; ++i) {
			fprintf(stdout, "%s%zd", (i == 0 ? "" : "\t"), result_a[i]);
		}
		fprintf(stdout, "\n");
		for (i = 0; i < result_b_count; ++i) {
			fprintf(stdout, "%s%zd", (i == 0 ? "" : "\t"), result_b[i]);
		}
		fprintf(stdout, "\n");
	}
	free(result_b);
	free(result_a);

	free(m);
	if (verbose) {
		printf("==============================\n");
	}
	return 0;
}

static int align_between_maps(const struct nick_map *map1, const struct nick_map *map2)
{
	size_t i, j, k;
	size_t map1_max_count = 0;
	size_t map2_max_count = 0;
	int *a, *b;

	for (i = 0; i < map1->size; ++i) {
		if (map1_max_count < map1->data[i].size) {
			map1_max_count = map1->data[i].size;
		}
	}
	++map1_max_count;

	for (i = 0; i < map2->size; ++i) {
		if (map2_max_count < map2->data[i].size) {
			map2_max_count = map2->data[i].size;
		}
	}
	++map2_max_count;

	a = malloc(sizeof(int) * map1_max_count);
	b = malloc(sizeof(int) * map2_max_count);
	if (!a || !b) {
		free(b);
		free(a);
		return -ENOMEM;
	}

	fprintf(stdout, "#>0\tAlignmentID\tMol0ID\tMol1ID\n");
	for (i = 0; i < map1->size; ++i) {
		const struct nick_list *list1 = &map1->data[i];
		a[0] = 0;
		for (k = 0; k < list1->size; ++k) {
			a[k + 1] = list1->data[k].pos;
		}

		for (j = (map1 == map2 ? i + 1: 0); j < map2->size; ++j) {
			const struct nick_list *list2 = &map2->data[j];
			b[0] = 0;
			for (k = 0; k < list2->size; ++k) {
				b[k + 1] = list2->data[k].pos;
			}

			align(map1->data[i].fragment_name, a, list1->size + 1,
				map2->data[j].fragment_name, b, list2->size + 1);
		}
	}
	free(b);
	free(a);
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
		ret = align_between_maps(&map, &map);
	} else {
		struct nick_map map2;
		nick_map_init(&map2);
		if (nick_map_load(&map2, argv[optind + 1])) {
			nick_map_free(&map2);
			nick_map_free(&map);
			return 1;
		}
		ret = align_between_maps(&map, &map2);
		nick_map_free(&map2);
	}
	nick_map_free(&map);
	return ret;
}
