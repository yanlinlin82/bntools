#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include "nick_map.h"
#include "version.h"

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

int nick_map_init(struct nick_map *map, const char *enzyme, const char *rec_seq)
{
	map->capacity = 0;
	map->size = 0;
	map->data = NULL;
	map->enzyme = strdup(enzyme);
	map->rec_seq = strdup(rec_seq);
	if (!map->enzyme || !map->rec_seq) {
		free(map->enzyme);
		free(map->rec_seq);
		return -ENOMEM;
	}
	return 0;
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

struct nick_list *nick_map_add_chrom(struct nick_map *map, const char *chrom)
{
	struct nick_list *p;
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
	memset(p, 0, sizeof(struct nick_list));
	p->chrom_name = strdup(chrom);
	if (!p->chrom_name) {
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

static void write_cmap_line(gzFile file, const char *cmap_id, int contig_length,
		size_t num_sites, size_t site_id, int label_channel, int position,
		int stddev, int coverage, int occurance)
{
	gzprintf(file, "%s\t%d\t%zd\t%zd\t%d\t%d\t%d\t%d\t%d\n",
			cmap_id, contig_length, num_sites, site_id,
			label_channel, position, stddev, coverage, occurance);
}

static void write_command_line(gzFile file)
{
	char name[64] = "";
	snprintf(name, sizeof(name), "/proc/%d/cmdline", getpid());
	FILE *fp = fopen(name, "r");
	if (fp) {
		gzprintf(file, "##commandline=");
		while (!feof(fp)) {
			int c = fgetc(fp);
			if (c == EOF) break;
			gzprintf(file, "%c", (c ? c : ' '));
		}
		gzprintf(file, "\n");
		fclose(fp);
	}
}

void nick_map_write(gzFile file, const struct nick_map *map, const struct nick_list *list)
{
	static const char * const STRAND[] = { "?", "+", "-", "+/-" };
	size_t i, j;
	for (i = 0; i < map->size; ++i) {
		const struct nick_list *p = &map->data[i];
		if (i == 0) {
			gzprintf(file, "##fileformat=MAPv0.1\n");
			gzprintf(file, "##enzyme=%s/%s\n", map->enzyme, map->rec_seq);
			gzprintf(file, "##program=bntools\n");
			gzprintf(file, "##programversion="VERSION"\n");
			write_command_line(file);
			gzprintf(file, "#chrom\tpos\tstrand\n");
		}
		if (!list || list == p) {
			for (j = 0; j < p->size; ++j) {
				const struct nick *q = &p->data[j];
				gzprintf(file, "%s\t%d\t%s\n",
						p->chrom_name, q->pos, STRAND[q->strand]);
			}
		}
	}
}

void nick_map_write_cmap(gzFile file, const struct nick_map *map)
{
	size_t i, j;

	gzprintf(file, "# CMAP File Version:  0.1\n");
	gzprintf(file, "# Label Channels:  1\n");
	gzprintf(file, "# Nickase Recognition Site 1:  %s/%s\n", map->enzyme, map->rec_seq);
	gzprintf(file, "# Number of Consensus Nanomaps:    %zd\n", map->size);
	gzprintf(file, "#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence\n");
	gzprintf(file, "#f int	float	int	int	int	float	float	int	int\n");

	for (i = 0; i < map->size; ++i) {
		const struct nick_list *p = &map->data[i];
		for (j = 0; j < p->size; ++j) {
			write_cmap_line(file, p->chrom_name, p->chrom_size,
					p->size, j, 1, p->data[j].pos, 0, 0, 0);
		}
		if (p->size > 0) {
			write_cmap_line(file, p->chrom_name, p->chrom_size,
					p->size, p->size + 1, 1, p->chrom_size, 0, 1, 1);
		}
	}
}
