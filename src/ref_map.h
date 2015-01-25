#ifndef __REF_MAP_H__
#define __REF_MAP_H__

#include <string.h>
#include "nick_map.h"

struct ref_map {
	struct nick_map map;

	char rec_bases[MAX_REC_SEQ_SIZE];
	int rec_seq_size;
	int nick_offset;
	int palindrome;
};

void ref_map_init(struct ref_map *ref);
void ref_map_free(struct ref_map *ref);
int ref_map_set_enzyme(struct ref_map *ref, const char *enzyme, const char *rec_seq);

int nick_map_load_fasta(struct ref_map *ref, const char *filename, int chrom_only, int verbose);

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

#endif /* __REF_MAP_H__ */
