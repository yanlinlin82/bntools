#ifndef __REF_MAP_H__
#define __REF_MAP_H__

#include <string.h>
#include "nick_map.h"
#include "array.h"

enum node_flag {
	FIRST_INTERVAL = 1,
	LAST_INTERVAL  = 2,
};

struct ref_node {
	size_t chrom;  /* item index in ref */
	size_t label;  /* label index in fragment/chrom */
	int pos;
	int size;
	int flag;
};

struct ref_index {
	const struct ref_node *node;
	int direct;
	int uniq_count;
};

struct ref_map {
	struct nick_map map;

	char rec_bases[MAX_REC_SEQ_SIZE];
	int rec_seq_size;
	int nick_offset;
	int palindrome;

	array(struct ref_node) nodes;
	array(struct ref_index) index_;
};

void ref_map_init(struct ref_map *ref);
void ref_map_free(struct ref_map *ref);
int ref_map_set_enzyme(struct ref_map *ref, const char *enzyme, const char *rec_seq);

int nick_map_load_fasta(struct ref_map *ref, const char *filename, int chrom_only, int verbose);

int ref_map_build_index(struct ref_map *ref);
void ref_map_dump(const struct ref_map *ref);

#endif /* __REF_MAP_H__ */
