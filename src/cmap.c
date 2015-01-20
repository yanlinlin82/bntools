#include <stdlib.h>
#include <errno.h>
#include "cmap.h"

#define MIN_INCREMENT 16
#define MAX_INCREMENT 128

size_t new_capacity(size_t capacity, size_t demand_size)
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

static int cmap_reserve(struct cmap *cmap, size_t size)
{
	if (cmap->capacity < size) {
		size_t capacity = new_capacity(cmap->capacity, size);
		struct nick_site_list *p = malloc(sizeof(struct nick_site_list) * capacity);
		if (!p) {
			return -ENOMEM;
		}
		if (cmap->size > 0) {
			memcpy(p, cmap->data, sizeof(struct nick_site_list) * cmap->size);
		}
		free(cmap->data);
		cmap->data = p;
		cmap->capacity = capacity;
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

struct nick_site_list *cmap_add_chrom(struct cmap *cmap, const char *chrom)
{
	struct nick_site_list *p;
	size_t i;

	for (i = 0; i < cmap->size; ++i) {
		if (strcmp(cmap->data[i].chrom_name, chrom) == 0) {
			return &cmap->data[i];
		}
	}
	if (cmap_reserve(cmap, i + 1)) {
		return NULL;
	}

	p = &cmap->data[i];
	memset(p, 0, sizeof(struct nick_site_list));
	p->chrom_name = strdup(chrom);
	if (!p->chrom_name) {
		return NULL;
	}
	++cmap->size;
	return p;
}

int cmap_add_site(struct nick_site_list *p, int pos, int strand)
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

void cmap_free(struct cmap *cmap)
{
	size_t i;
	for (i = 0; i < cmap->size; ++i) {
		free(cmap->data[i].chrom_name);
		free(cmap->data[i].data);
	}
	free(cmap->data);
	cmap->data = NULL;
	cmap->size = 0;
	cmap->capacity= 0;
}
