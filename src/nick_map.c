#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include "nick_map.h"
#include "base_map.h"

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
		struct nick_list *p = malloc(sizeof(struct nick_list) * capacity);
		if (!p) {
			return -ENOMEM;
		}
		if (map->size > 0) {
			memcpy(p, map->data, sizeof(struct nick_list) * map->size);
		}
		free(map->data);
		map->data = p;
		map->capacity = capacity;
	}
	return 0;
}

static int nick_list_reserve(struct nick_list *list, size_t size)
{
	if (list->capacity < size) {
		size_t capacity = new_capacity(list->capacity, size);
		struct nick *p = malloc(sizeof(struct nick) * capacity);
		if (!p) {
			return -ENOMEM;
		}
		if (list->size > 0) {
			memcpy(p, list->data, sizeof(struct nick) * list->size);
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
	map->nick_offset = -1;
}

void nick_map_free(struct nick_map *map)
{
	size_t i;

	for (i = 0; i < map->size; ++i) {
		free(map->data[i].fragment_name);
		free(map->data[i].data);
	}
	free(map->data);
	map->data = NULL;
	map->size = 0;
	map->capacity= 0;
}

int nick_map_set_enzyme(struct nick_map *map, const char *enzyme, const char *rec_seq)
{
	char checked_rec_seq[sizeof(map->rec_seq)];
	char rec_bases[sizeof(map->rec_bases)];
	int rec_seq_size;
	int nick_offset = -1;
	int palindrome;
	int i;

	for (i = 0, rec_seq_size = 0; rec_seq[i]; ++i) {
		if (rec_seq[i] == '^') {
			if (nick_offset >= 0) {
				fprintf(stderr, "Error: Invalid recognition sequence '%s'\n", rec_seq);
				return 1;
			}
			nick_offset = i;
		} else {
			char c = char_to_base(rec_seq[i]);
			if (c == 0) {
				fprintf(stderr, "Error: Invalid character '%c' in recognition sequence '%s'\n", rec_seq[i], rec_seq);
				return 1;
			}
			if (rec_seq_size >= sizeof(rec_bases)) {
				fprintf(stderr, "Error: Recognition sequence is too long\n");
				return 1;
			}
			rec_bases[rec_seq_size++] = c;
		}
	}
	for (i = 0; i < rec_seq_size / 2; ++i) {
		if (rec_bases[i] != base_to_comp(rec_bases[rec_seq_size - i - 1])) {
			palindrome = 0;
			break;
		}
	}
	snprintf(checked_rec_seq, sizeof(checked_rec_seq), "%s%s",
			(nick_offset < 0 ? "^" : ""), rec_seq);

	if (map->nick_offset >= 0) { /* merging maps */
		if (strcmp(map->enzyme, enzyme) != 0 ||
				strcmp(map->rec_seq, checked_rec_seq) != 0) {
			fprintf(stderr, "Error: Merging is supported for only single enzyme\n");
			return 1;
		}
	} else {
		snprintf(map->enzyme, sizeof(map->enzyme), "%s", enzyme);
		snprintf(map->rec_seq, sizeof(map->rec_seq), "%s", checked_rec_seq);
		assert(strcmp(map->enzyme, enzyme) == 0);
		assert(strcmp(map->rec_seq, checked_rec_seq) == 0);

		memcpy(map->rec_bases, rec_bases, sizeof(map->rec_bases));
		map->rec_seq_size = rec_seq_size;
		map->nick_offset = (nick_offset < 0 ? 0 : nick_offset);
		map->palindrome = palindrome;
	}
	return 0;
}

struct nick_list *nick_map_add_fragment(struct nick_map *map, const char *name)
{
	struct nick_list *p;
	size_t i;

	for (i = 0; i < map->size; ++i) {
		if (strcmp(map->data[i].fragment_name, name) == 0) {
			return &map->data[i];
		}
	}
	if (nick_map_reserve(map, i + 1)) {
		return NULL;
	}

	p = &map->data[i];
	memset(p, 0, sizeof(struct nick_list));
	p->fragment_name = strdup(name);
	if (!p->fragment_name) {
		return NULL;
	}
	++map->size;
	return p;
}

int nick_map_add_site(struct nick_list *p, int pos, int strand)
{
	size_t i, j;

	for (i = p->size; i > 0; --i) {
		if (p->data[i - 1].pos == pos) {
			p->data[i - 1].strand |= strand;
			return 0;
		} else if (p->data[i - 1].pos < pos) {
			break;
		}
	}
	if (nick_list_reserve(p, p->size + 1)) {
		return -ENOMEM;
	}
	for (j = p->size; j > i; --j) {
		memcpy(&p->data[j], &p->data[j - 1], sizeof(struct nick));
	}
	p->data[i].pos = pos;
	p->data[i].strand = strand;
	++p->size;
	return 0;
}
