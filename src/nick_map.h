#ifndef __NICK_MAP_H__
#define __NICK_MAP_H__

#include <string.h>
#include <zlib.h>

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
	char *chrom_name;
	int chrom_size;
};

struct nick_map {
	size_t capacity;
	size_t size;
	struct nick_list *data;
	char *enzyme;
	char *rec_seq;
};

void nick_map_init(struct nick_map *map);
void nick_map_free(struct nick_map *map);

int nick_map_set_enzyme(struct nick_map *map, const char *enzyme, const char *rec_seq);

struct nick_list *nick_map_add_chrom(struct nick_map *map, const char *chrom);
int nick_map_add_site(struct nick_list *p, int pos, int strand);

int nick_map_load(struct nick_map *map, gzFile file);

void nick_map_write(gzFile file, const struct nick_map *map, const struct nick_list *list);
void nick_map_write_cmap(gzFile file, const struct nick_map *map);

#endif /* __NICK_MAP_H__ */
