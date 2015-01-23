#ifndef __NICK_MAP_H__
#define __NICK_MAP_H__

#include <string.h>

#define MAX_ENZYME_NAME_SIZE 16
#define MAX_REC_SEQ_SIZE 32
#define MAX_CHROM_NAME_SIZE 64

#define STRAND_UNKNOWN 0
#define STRAND_PLUS    1
#define STRAND_MINUS   2
#define STRAND_BOTH    (STRAND_PLUS | STRAND_MINUS)
#define STRAND_END     4

struct nick {
	int pos;
	int strand;
};

struct nick_list {
	size_t capacity;
	size_t size;
	struct nick *data;
	char *fragment_name;
	int fragment_size;
};

struct nick_map {
	size_t capacity;
	size_t size;
	struct nick_list *data;

	char enzyme[MAX_ENZYME_NAME_SIZE];
	char rec_seq[MAX_REC_SEQ_SIZE];

	char rec_bases[MAX_REC_SEQ_SIZE];
	int rec_seq_size;
	int nick_offset;
	int palindrome;
};

void nick_map_init(struct nick_map *map);
void nick_map_free(struct nick_map *map);

int nick_map_set_enzyme(struct nick_map *map, const char *enzyme, const char *rec_seq);

struct nick_list *nick_map_add_fragment(struct nick_map *map, const char *name);
int nick_map_add_site(struct nick_list *p, int pos, int strand);

int nick_map_load(struct nick_map *map, const char *filename);
int nick_map_load_fasta(struct nick_map *map, const char *filename,
		int only_chromosome, int transform_to_number, int verbose);
int nick_map_save(const struct nick_map *map, const char *filename, int output_cmap);

#endif /* __NICK_MAP_H__ */
