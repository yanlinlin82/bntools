#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include "nick_map.h"
#include "base_map.h"

void nick_map_init(struct nick_map *map)
{
	memset(map, 0, sizeof(struct nick_map));
}

void nick_map_free(struct nick_map *map)
{
	size_t i;
	for (i = 0; i < map->fragments.size; ++i) {
		array_free(map->fragments.data[i].nicks);
	}
	array_free(map->fragments);
}

void nick_map_set_enzyme(struct nick_map *map, const char *enzyme, const char *rec_seq)
{
	snprintf(map->enzyme, sizeof(map->enzyme), "%s", enzyme);
	snprintf(map->rec_seq, sizeof(map->rec_seq), "%s", rec_seq);
	assert(strcmp(map->enzyme, enzyme) == 0);
	assert(strcmp(map->rec_seq, rec_seq) == 0);
}

struct fragment *nick_map_add_fragment(struct nick_map *map, const char *name)
{
	struct fragment *f;
	size_t i;

	for (i = 0; i < map->fragments.size; ++i) {
		if (strncmp(map->fragments.data[i].name, name, MAX_FRAGMENT_NAME_SIZE) == 0) {
			return &map->fragments.data[i];
		}
	}
	if (array_reserve(map->fragments, map->fragments.size + 1)) {
		return NULL;
	}

	f = &map->fragments.data[map->fragments.size++];
	memset(f, 0, sizeof(struct fragment));
	snprintf(f->name, sizeof(f->name), "%s", name);
	return f;
}

int nick_map_add_site(struct fragment *f, int pos, unsigned int flag)
{
	size_t i, j;

	for (i = f->nicks.size; i > 0; --i) {
		if (f->nicks.data[i - 1].pos == pos) {
			f->nicks.data[i - 1].flag |= flag;
			return 0;
		} else if (f->nicks.data[i - 1].pos < pos) {
			break;
		}
	}
	if (array_reserve(f->nicks, f->nicks.size + 1)) {
		return -ENOMEM;
	}
	for (j = f->nicks.size; j > i; --j) {
		memcpy(&f->nicks.data[j], &f->nicks.data[j - 1], sizeof(struct nick));
	}
	f->nicks.data[i].pos = pos;
	f->nicks.data[i].flag = flag;
	++f->nicks.size;
	return 0;
}
