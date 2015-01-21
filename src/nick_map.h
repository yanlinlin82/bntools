#ifndef __NICK_MAP_H__
#define __NICK_MAP_H__

#include <string.h>

#define STRAND_PLUS 0
#define STRAND_MINUS 1

struct nick_site {
	int pos;
	int strand;
};

struct nick_site_list {
	size_t capacity;
	size_t size;
	struct nick_site *data;
	char *chrom_name;
	int chrom_size;
};

struct nick_map {
	size_t capacity;
	size_t size;
	struct nick_site_list *data;
};

void nick_map_init(struct nick_map *map);
void nick_map_free(struct nick_map *map);

struct nick_site_list *nick_map_add_chrom(struct nick_map *nick_map, const char *chrom);
int nick_map_add_site(struct nick_site_list *p, int pos, int strand);

#endif /* __NICK_MAP_H__ */
