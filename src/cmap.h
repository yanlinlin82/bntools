#ifndef __CMAP_H__
#define __CMAP_H__

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

struct cmap {
	size_t capacity;
	size_t size;
	struct nick_site_list *data;
};

struct nick_site_list *cmap_add_chrom(struct cmap *cmap, const char *chrom);
int cmap_add_site(struct nick_site_list *p, int pos, int strand);
void cmap_free(struct cmap *cmap);

#endif /* __CMAP_H__ */
