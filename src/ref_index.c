#include <stdio.h>
#include <stdlib.h>
#include "ref_index.h"

size_t ref_node_count = 0;
struct node * ref_nodes = NULL;
struct node **ref_index = NULL;

static int sort_by_size(const void *a, const void *b)
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

void generate_ref_nodes(const struct nick_map *ref)
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
			ref_nodes[k].uniq_count = 0;
			ref_index[k] = &ref_nodes[k];
			++k;
		}
	}

	qsort(ref_index, ref_node_count, sizeof(struct node *), sort_by_size);

	for (i = 0; i + 1 < ref_node_count; ++i) {
		struct node *a = ref_index[i];
		struct node *b = ref_index[i + 1];
		struct node *end = ref_nodes + ref_node_count;
		for (j = 0; a + j < end && b + j < end; ++j) {
			if (a[j].size != b[j].size) {
				break;
			}
		}
		if (a->uniq_count < j + 1) {
			a->uniq_count = j + 1;
		}
		if (b->uniq_count < j + 1) {
			b->uniq_count = j + 1;
		}
	}
}

void print_sorted_ref(const struct nick_map *ref)
{
	size_t i, j;

	fprintf(stdout, "#id\tname\tpos\tsize\tflag\tuniq\tseq\n");

	for (i = 0; i < ref_node_count; ++i) {
		fprintf(stdout, "%zd\t%s\t%d\t%d\t%d\t%d\t", i + 1,
				ref->data[ref_index[i]->chrom].fragment_name,
				ref_index[i]->pos, ref_index[i]->size, ref_index[i]->flag,
				ref_index[i]->uniq_count);
		for (j = 0; j < ref_index[i]->uniq_count; ++j) {
			fprintf(stdout, "%s%d", (j == 0 ? "": ","), ref_index[i][j].size);
		}
		fprintf(stdout, "\n");
	}
}
