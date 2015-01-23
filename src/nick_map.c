#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include "nick_map.h"
#include "base_map.h"

void nick_map_init(struct nick_map *map)
{
	memset(map, 0, sizeof(struct nick_map));
	map->nick_offset = -1;
}

void nick_map_free(struct nick_map *map)
{
	size_t i;
	for (i = 0; i < map->fragments.size; ++i) {
		free(map->fragments.data[i].name);
		array_free(map->fragments.data[i].nicks);
	}
	array_free(map->fragments);
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

struct fragment *nick_map_add_fragment(struct nick_map *map, const char *name)
{
	struct fragment *f;
	size_t i;

	for (i = 0; i < map->fragments.size; ++i) {
		if (strcmp(map->fragments.data[i].name, name) == 0) {
			return &map->fragments.data[i];
		}
	}
	if (array_reserve(map->fragments, i + 1)) {
		return NULL;
	}

	f = &map->fragments.data[i];
	memset(f, 0, sizeof(struct fragment));
	f->name = strdup(name);
	if (!f->name) {
		return NULL;
	}
	++map->fragments.size;
	return f;
}

int nick_map_add_site(struct fragment *f, int pos, int strand)
{
	size_t i, j;

	for (i = f->nicks.size; i > 0; --i) {
		if (f->nicks.data[i - 1].pos == pos) {
			f->nicks.data[i - 1].strand |= strand;
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
	f->nicks.data[i].strand = strand;
	++f->nicks.size;
	return 0;
}
