#ifndef __REF_INDEX_H__
#define __REF_INDEX_H__

#include <string.h>
#include "nick_map.h"

enum node_flag {
	FIRST_INTERVAL = 1,
	LAST_INTERVAL  = 2,
};

struct node {
	size_t chrom;  /* item index in ref */
	int pos;
	int size;
	int flag;
	int uniq_count;
};

extern size_t ref_node_count;
extern struct node * ref_nodes;
extern struct node **ref_index;

void generate_ref_nodes(const struct nick_map *ref);
void print_sorted_ref(const struct nick_map *ref);

#endif /* __REF_INDEX_H__ */
