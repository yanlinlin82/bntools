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

	array(struct ref_node) nodes;
	array(struct ref_index) index_;
};

void ref_map_init(struct ref_map *ref);
void ref_map_free(struct ref_map *ref);

struct rec_site {
	char enzyme[MAX_ENZYME_NAME_SIZE + 1];
	char rec_seq[MAX_REC_SEQ_SIZE + 1];
	char rec_bases[MAX_REC_SEQ_SIZE];
	int rec_seq_size;
	int nick_offset;
	int palindrome;
};

int prepare_rec_site(struct rec_site *site, const char *enzyme, const char *rec_seq);

int nick_map_load_seq(struct ref_map *ref, const char *filename,
		const struct rec_site *site, int chrom_only, int verbose);

int ref_map_build_index(struct ref_map *ref);
int ref_map_save(const struct ref_map *ref, const char *filename);
int ref_map_load(struct ref_map *ref, const char *filename);

const char *get_index_filename(const char *filename, char *buf, size_t bufsize);

#endif /* __REF_MAP_H__ */
