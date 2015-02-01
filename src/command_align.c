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

static int verbose = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools align [options] <map_a> [<map_b>]\n"
			"\n"
			"Options:\n"
			"   <map_a/b>   input map file(s), in tsv/cmap/bnx format\n"
			"   -v          show verbose message\n"
			"   -h          show this help\n"
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

static int align(const struct fragment *fa, const struct fragment *fb)
{
	struct node *m;
	size_t size, i, j, delta_i, delta_j;
	size_t *result_a, *result_b, result_a_count, result_b_count;
	size_t w = fa->nicks.size;
	size_t h = fb->nicks.size;

	size = sizeof(struct node) * w * h;
	m = malloc(size);
	if (!m) {
		fprintf(stderr, "Error: Failed to allocate memory for DP matrix!\n");
		return -ENOMEM;
	}
	memset(m, 0, size);

	if (verbose) {
		fprintf(stdout, "align between '%s' and '%s'\n", fa->name, fb->name);
	}

	for (j = 1; j < h; ++j) {
		for (i = 1; i < w; ++i) {
			size_t index = j * w + i;

			for (delta_j = 1; delta_j < j; ++j) {
				int frag_b = fb->nicks.data[j].pos - fb->nicks.data[j - delta_j].pos;

				for (delta_i = 1; delta_i < i; ++i) {
					int frag_a = fa->nicks.data[i].pos - fa->nicks.data[i - delta_i].pos;

					if (delta_i > 1 || delta_j > 1) continue;
					if (frag_a > frag_b * (1 + TOLERANCE)) break;

					if (abs(frag_a - frag_b) < min(frag_a, frag_b) * TOLERANCE) {
						int score = m[(j - delta_j) * w + (i - delta_i)].score + 1;
						if (m[index].score < score) {
							m[index].score = score;
							m[index].delta_i = delta_i;
							m[index].delta_j = delta_j;
						}
					}
				}
			}
		}
	}

	if (verbose > 0) {
		printf("a (%zd):\n", w);
		for (i = 1; i < w; ++i) {
			printf(" [%2zd] %-5d", i, fa->nicks.data[i].pos - fa->nicks.data[i - 1].pos);
			if (i % 5 == 0) printf("\n");
		}
		if (i % 5 != 1) printf("\n");

		printf("b (%zd):\n", h);
		for (i = 1; i < h; ++i) {
			printf(" [%2zd] %-5d", i, fb->nicks.data[i].pos - fb->nicks.data[i - 1].pos);
			if (i % 5 == 0) printf("\n");
		}
		if (i % 5 != 1) printf("\n");

		printf("matrix:\n");
		for (j = 1; j < h; ++j) {
			for (i = 1; i < w; ++i) {
				const struct node *p = &m[j * w + i];
				printf("%2d/%d,%d", p->score, p->delta_i, p->delta_j);
			}
			printf("\n");
		}
	}

	result_a = malloc(sizeof(size_t) * w);
	result_b = malloc(sizeof(size_t) * h);
	for (;;) {
		int max_score = 0;
		size_t max_score_i = 0;
		size_t max_score_j = 0;

		for (j = 0; j < h; ++j) {
			for (i = 0; i < w; ++i) {
				if (max_score < m[j * w + i].score) {
					max_score = m[j * w + i].score;
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
			struct node *p = &m[j * w + i];

			result_a[result_a_count++] = i;
			result_b[result_b_count++] = j;

			if (p->delta_i == 0 && p->delta_j == 0) break;

			if (verbose > 0) {
				printf("a: %d { ", fa->nicks.data[i].pos - fa->nicks.data[i - p->delta_i].pos);
				for (k = 0; k < p->delta_i; ++k) {
					size_t index = i - p->delta_i + k;
					if (k > 0) printf(", ");
					printf("[%zd] %d", index, fa->nicks.data[index + 1].pos - fa->nicks.data[index].pos);
				}
				printf(" }\t");

				printf("b: %d { ", fb->nicks.data[i].pos - fb->nicks.data[i - p->delta_j].pos);
				for (k = 0; k < p->delta_j; ++k) {
					size_t index = j - p->delta_j + k;
					if (k > 0) printf(", ");
					printf("[%zd] %d", index, fb->nicks.data[index + 1].pos - fb->nicks.data[index].pos);
				}
				printf(" }\n");
			}

			i -= p->delta_i;
			j -= p->delta_j;
			p->score = 0;
		}

		{
			static int count = 0;
			fprintf(stdout, ">0\t%d\t%s\t%s\n", ++count, fa->name, fb->name);
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
	size_t i, j;
	fprintf(stdout, "#>0\tAlignmentID\tMol0ID\tMol1ID\n");
	for (i = 0; i < map1->fragments.size; ++i) {
		for (j = (map1 == map2 ? i + 1: 0); j < map2->fragments.size; ++j) {
			align(&map1->fragments.data[i], &map2->fragments.data[j]);
		}
	}
	return 0;
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
