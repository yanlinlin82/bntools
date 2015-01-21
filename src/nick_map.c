#include <stdlib.h>
#include <errno.h>
#include "nick_map.h"

#define MIN_INCREMENT 16
#define MAX_INCREMENT 128

static size_t new_capacity(size_t capacity, size_t demand_size)
{
	while (capacity < demand_size) {
		if (capacity < MIN_INCREMENT) {
			capacity += MIN_INCREMENT;
		} else if (capacity> MAX_INCREMENT) {
			capacity += MAX_INCREMENT;
		} else {
			capacity *= 2;
		}
	}
	return capacity;
}

static int nick_map_reserve(struct nick_map *map, size_t size)
{
	if (map->capacity < size) {
		size_t capacity = new_capacity(map->capacity, size);
		struct nick_site_list *p = malloc(sizeof(struct nick_site_list) * capacity);
		if (!p) {
			return -ENOMEM;
		}
		if (map->size > 0) {
			memcpy(p, map->data, sizeof(struct nick_site_list) * map->size);
		}
		free(map->data);
		map->data = p;
		map->capacity = capacity;
	}
	return 0;
}

static int nick_site_list_reserve(struct nick_site_list *list, size_t size)
{
	if (list->capacity < size) {
		size_t capacity = new_capacity(list->capacity, size);
		struct nick_site *p = malloc(sizeof(struct nick_site) * capacity);
		if (!p) {
			return -ENOMEM;
		}
		if (list->size > 0) {
			memcpy(p, list->data, sizeof(struct nick_site) * list->size);
		}
		free(list->data);
		list->data = p;
		list->capacity = capacity;
	}
	return 0;
}

void nick_map_init(struct nick_map *map)
{
	memset(map, 0, sizeof(struct nick_map));
}

void nick_map_free(struct nick_map *map)
{
	size_t i;
	for (i = 0; i < map->size; ++i) {
		free(map->data[i].chrom_name);
		free(map->data[i].data);
	}
	free(map->data);
	map->data = NULL;
	map->size = 0;
	map->capacity= 0;
}

struct nick_site_list *nick_map_add_chrom(struct nick_map *map, const char *chrom)
{
	struct nick_site_list *p;
	size_t i;

	for (i = 0; i < map->size; ++i) {
		if (strcmp(map->data[i].chrom_name, chrom) == 0) {
			return &map->data[i];
		}
	}
	if (nick_map_reserve(map, i + 1)) {
		return NULL;
	}

	p = &map->data[i];
	memset(p, 0, sizeof(struct nick_site_list));
	p->chrom_name = strdup(chrom);
	if (!p->chrom_name) {
		return NULL;
	}
	++map->size;
	return p;
}

int nick_map_add_site(struct nick_site_list *p, int pos, int strand)
{
	size_t i;

	if (nick_site_list_reserve(p, p->size + 1)) {
		return -ENOMEM;
	}

	for (i = p->size; i > 0; --i) {
		if (p->data[i - 1].pos > pos) {
			memcpy(&p->data[i], &p->data[i - 1], sizeof(struct nick_site));
		} else {
			break;
		}
	}
	p->data[i].pos = pos;
	p->data[i].strand = strand;
	++p->size;
	return 0;
}
