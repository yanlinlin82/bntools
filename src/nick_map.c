#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <zlib.h>
#include "base_map.h"
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

void nick_map_init(struct nick_map *map)
{
	memset(map, 0, sizeof(struct nick_map));
	map->nick_offset = -1;
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

static inline int string_begins_as(const char *s, const char *prefix)
{
	return (memcmp(s, prefix, strlen(prefix)) == 0);
}

static inline int skip_to_next_line(gzFile file, char *buf, size_t bufsize)
{
	while (strchr(buf, '\n') == NULL) {
		if (!gzgets(file, buf, bufsize)) return 1;
	}
	return 0;
}

static int seq_match(const char *ref, const char *query, size_t len, int strand)
{
	size_t i;
	assert(strand == STRAND_PLUS || strand == STRAND_MINUS);
	for (i = 0; i < len; ++i) {
		char r = ref[i];
		char q = (strand == STRAND_PLUS ? query[i] : base_to_comp(query[len - i - 1]));
		if ((r & q) != r) {
			return 0;
		}
	}
	return 1;
}

struct buffer {
	char data[MAX_REC_SEQ_SIZE * 2];
	size_t pos;
};

static int process_line(struct nick_map *map, struct nick_list *list,
		const char *line, const char *chrom, int base_count, struct buffer *buf)
{
	const char *p;
	int strand;
	int matched;
	for (p = line; *p; ++p) {
		if (isspace(*p)) {
			continue;
		}
		if (buf->pos >= sizeof(buf->data)) {
			memcpy(buf->data, buf->data + sizeof(buf->data) - map->rec_seq_size + 1, map->rec_seq_size - 1);
			buf->pos = map->rec_seq_size - 1;
		}
		++base_count;
		buf->data[buf->pos++] = char_to_base(*p);
		if (buf->pos < map->rec_seq_size) {
			continue;
		}
		for (strand = STRAND_PLUS, matched = 0; strand <= STRAND_MINUS; ++strand) {
			if (matched || seq_match(buf->data + buf->pos - map->rec_seq_size, map->rec_bases, map->rec_seq_size, strand)) {
				int site_pos = base_count - (strand == STRAND_MINUS ? map->nick_offset : (map->rec_seq_size - map->nick_offset));
				if (nick_map_add_site(list, site_pos, strand)) {
					return -ENOMEM;
				}
				matched = map->palindrome;
			}
		}
	}
	return base_count;
}

static int process_map(gzFile fin, struct nick_map *map, int transform_to_number, int verbose)
{
	struct nick_list *list = NULL;
	char chrom[MAX_CHROM_NAME_SIZE] = "";
	int base_count = 0;
	struct buffer buf = { };

	while (!gzeof(fin)) {
		char line[256];
		if (!gzgets(fin, line, sizeof(line))) break;
		if (line[0] == '>') {
			if (list) {
				list->chrom_size = base_count;
				if (verbose > 0) {
					fprintf(stderr, "%d bp\n", base_count);
				}
			}
			if (transform_to_number) {
				static int number = 0;
				snprintf(chrom, sizeof(chrom), "%d", ++number);
			} else {
				char *p = line + 1;
				while (*p && !isspace(*p)) ++p;
				*p = '\0';
				snprintf(chrom, sizeof(chrom), line + 1);
			}
			if (verbose > 0) {
				fprintf(stderr, "Loading sequence '%s' ... ", chrom);
			}
			base_count = 0;
			list = nick_map_add_chrom(map, chrom);
			if (!list) {
				return -ENOMEM;
			}
		} else {
			int n = process_line(map, list, line, chrom, base_count, &buf);
			if (n < 0) {
				return n;
			}
			base_count = n;
		}
	}
	if (list) {
		list->chrom_size = base_count;
		if (verbose > 0) {
			fprintf(stderr, "%d bp\n", base_count);
		}
	}
	return 0;
}

int nick_map_load_fasta(struct nick_map *map, const char *filename, int transform_to_number, int verbose)
{
	gzFile file;
	int c;

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		file = gzdopen(0, "r"); /* stdin */
	} else {
		file = gzopen(filename, "r");
	}
	if (!file) {
		fprintf(stderr, "Error: Can not open FASTA file '%s'\n", filename);
		return 1;
	}

	c = gzgetc(file);
	if (c != '>') {
		fprintf(stderr, "Error: File '%s' is not in FASTA format\n", filename);
		gzclose(file);
		return 1;
	}
	gzungetc(c, file);

	if (process_map(file, map, transform_to_number, verbose)) {
		return 1;
	}

	gzclose(file);
	return 0;
}

static inline int to_integer(double x) { return (int)(x + .5); }

static int load_bnx(gzFile file, long long lineNo, struct nick_map *map)
{
	char buf[256];
	char molecule_id[64] = "";
	double length = 0;
	int molecule_length;
	struct nick_list *list = NULL;
	while (!gzeof(file)) {
		++lineNo;
		if (!gzgets(file, buf, sizeof(buf))) break;
		if (buf[0] == '#') continue;
		if (memcmp(buf, "0\t", 2) == 0) { /* molecule info */
			if (sscanf(buf + 2, "%s%lf", molecule_id, &length) != 2) {
				fprintf(stderr, "Error: Invalid format on line %lld\n", lineNo);
				return 1;
			}
			molecule_length = to_integer(length);
			list = nick_map_add_chrom(map, molecule_id);
			list->chrom_size = molecule_length;
			if (skip_to_next_line(file, buf, sizeof(buf))) break;
		} else if (memcmp(buf, "1\t", 2) == 0) { /* label positions */
			char *p = buf + 2;
			for (;;) {
				char *q = p;
				double value;
				while (*q && *q != '\t' && *q != '\n') ++q;
				if (*q == '\0') {
					size_t size = q - p;
					memcpy(buf, p, size);
					if (!gzgets(file, buf + size, sizeof(buf) - size)) break;
					p = buf + size;
					continue;
				}
				if (sscanf(p, "%lf", &value) != 1) {
					fprintf(stderr, "Error: Failed in reading float value on line %lld\n", lineNo);
					return 1;
				}
				nick_map_add_site(list, to_integer(value), 0);
				if (*q == '\t') {
					p = q + 1;
				} else {
					assert(*q == '\n');
					break;
				}
			}
		} else {
			if (skip_to_next_line(file, buf, sizeof(buf))) break;
		}
	}
	return 0;
}

static int load_cmap(gzFile file, long long lineNo, struct nick_map *map)
{
	char buf[256];
	struct nick_list *list = NULL;
	char lastMapId[32] = "";

	while (!gzeof(file)) {
		++lineNo;
		if (!gzgets(file, buf, sizeof(buf))) break;
		if (buf[0] == '#') {
			if (string_begins_as(buf, "# Nickase Recognition Site 1:")) {
				char *p, *q, *e;
				p = strchr(buf, ':');
				assert(p != NULL);
				++p;
				while (*p && isblank(*p)) ++p;
				q = strchr(p, '/');
				if (q != NULL) {
					*q++ = '\0';
					e = strchr(q, '\n');
					if (e) {
						*e = '\0';
					}
					if (nick_map_set_enzyme(map, p, q)) {
						return 1;
					}
				}
			}
			continue;
		} else {
			char mapId[32];
			int ctgLen;
			int numSites;
			int siteId;
			int labelChannel;
			int position;
			if (sscanf(buf, "%s%d%d%d%d%d", mapId, &ctgLen, &numSites, &siteId, &labelChannel, &position) != 6) {
				fprintf(stderr, "Error: Failed to parse data on line %lld\n", lineNo);
				return -EINVAL;
			}

			if (strcmp(lastMapId, mapId) != 0) {
				list = nick_map_add_chrom(map, mapId);
			}
			assert(list != NULL);

			if (labelChannel == 1) {
				nick_map_add_site(list, position, 0);
			} else {
				assert(labelChannel == 0);
				list->chrom_size = position;
			}
		}
	}
	return 0;
}

static int load_tsv(gzFile file, long long lineNo, struct nick_map *map)
{
	char buf[256];
	struct nick_list *list = NULL;
	char lastChrom[32] = "";

	while (!gzeof(file)) {
		++lineNo;
		if (!gzgets(file, buf, sizeof(buf))) break;
		if (buf[0] == '#') {
			if (string_begins_as(buf, "##enzyme=")) {
				char *p, *q, *e;
				p = strchr(buf, '=');
				assert(p != NULL);
				++p;
				q = strchr(p, '/');
				if (q != NULL) {
					*q++ = '\0';
					e = strchr(q, '\n');
					if (e) {
						*e = '\0';
					}
					if (nick_map_set_enzyme(map, p, q)) {
						return 1;
					}
				}
			}
			continue;
		} else {
			char chrom[32];
			int pos;
			char strandText[8];
			int strand;

			if (sscanf(buf, "%s%d%s", chrom, &pos, strandText) != 3) {
				fprintf(stderr, "Error: Failed to parse data on line %lld\n", lineNo);
				return -EINVAL;
			}

			if (strcmp(strandText, "?") == 0) {
				strand = STRAND_UNKNOWN;
			} else if (strcmp(strandText, "+") == 0) {
				strand = STRAND_PLUS;
			} else if (strcmp(strandText, "-") == 0) {
				strand = STRAND_MINUS;
			} else if (strcmp(strandText, "+/-") == 0) {
				strand = STRAND_BOTH;
			} else if (strcmp(strandText, "*") == 0) {
				strand = STRAND_END;
			} else {
				fprintf(stderr, "Error: Unknown strand text '%s' on line %lld\n", strandText, lineNo);
				return -EINVAL;
			}

			if (strcmp(lastChrom, chrom) != 0) {
				list = nick_map_add_chrom(map, chrom);
			}
			assert(list != NULL);

			if (strand != STRAND_END) {
				nick_map_add_site(list, pos, strand);
			} else {
				list->chrom_size = pos;
			}
		}
	}
	return 0;
}

int nick_map_load(struct nick_map *map, const char *filename)
{
	long long lineNo = 0;
	char buf[256];
	int format = 0; /* 1. BNX; 2. CMAP; 3. TSV */
	gzFile file;
	int ret = 0;

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		file = gzdopen(0, "r"); /* stdin */
	} else {
		file = gzopen(filename, "r");
	}
	if (!file) {
		fprintf(stderr, "Error: Can not open FASTA file '%s'\n", filename);
		return 1;
	}

	while (!gzeof(file)) {
		++lineNo;
		if (!gzgets(file, buf, sizeof(buf))) break;
		if (buf[0] != '#') break;
		if (string_begins_as(buf, "# BNX File Version:")) {
			format = 1;
			break;
		} else if (string_begins_as(buf, "# CMAP File Version:")) {
			format = 2;
			break;
		} else if (string_begins_as(buf, "##fileformat=MAPv0.1")) {
			format = 3;
			break;
		}
		if (skip_to_next_line(file, buf, sizeof(buf))) break;
	}
	if (format == 1) {
		ret = load_bnx(file, lineNo, map);
	} else if (format == 2) {
		ret = load_cmap(file, lineNo, map);
	} else if (format == 3) {
		ret = load_tsv(file, lineNo, map);
	} else {
		fprintf(stderr, "Error: Unknown input map format!\n");
		ret = 1;
	}
	gzclose(file);
	return ret;
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

static void nick_map_write(gzFile file, const struct nick_map *map)
{
	static const char * const STRAND[] = { "?", "+", "-", "+/-", "*" };
	size_t i, j;

	gzprintf(file, "##fileformat=MAPv0.1\n");
	gzprintf(file, "##enzyme=%s/%s\n", map->enzyme, map->rec_seq);
	gzprintf(file, "##program=bntools\n");
	gzprintf(file, "##programversion="VERSION"\n");
	write_command_line(file);
	gzprintf(file, "#chrom\tpos\tstrand\tsize\n");

	for (i = 0; i < map->size; ++i) {
		const struct nick_list *p = &map->data[i];
		for (j = 0; j < p->size; ++j) {
			const struct nick *q = &p->data[j];
			gzprintf(file, "%s\t%d\t%s\t%d\n",
					p->chrom_name, q->pos, STRAND[q->strand],
					q->pos - (j == 0 ? 0 : p->data[j - 1].pos));
		}
		gzprintf(file, "%s\t%d\t*\t%d\n", p->chrom_name, p->chrom_size,
				p->chrom_size - (p->size == 0 ? 0 : p->data[p->size - 1].pos));
	}
}

static void nick_map_write_cmap(gzFile file, const struct nick_map *map)
{
	size_t i, j;

	gzprintf(file, "# CMAP File Version:  0.1\n");
	gzprintf(file, "# Label Channels:  1\n");
	gzprintf(file, "# Nickase Recognition Site 1:  %s/%s\n", map->enzyme, map->rec_seq);
	gzprintf(file, "# Number of Consensus Nanomaps:    %zd\n", map->size);
	gzprintf(file, "#h CMapId\tContigLength\tNumSites\tSiteID"
			"\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n");
	gzprintf(file, "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n");

	for (i = 0; i < map->size; ++i) {
		const struct nick_list *p = &map->data[i];
		for (j = 0; j < p->size; ++j) {
			write_cmap_line(file, p->chrom_name, p->chrom_size,
					p->size, j, 1, p->data[j].pos, 0, 0, 0);
		}
		write_cmap_line(file, p->chrom_name, p->chrom_size,
				p->size, p->size + 1, 0, p->chrom_size, 0, 1, 1);
	}
}

int nick_map_save(const struct nick_map *map, const char *filename, int output_cmap)
{
	gzFile file;

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdout") == 0) {
		file = gzdopen(1, "wT"); /* stdout, without compression */
	} else {
		size_t len = strlen(filename);
		if (len > 3 && strcmp(filename + len - 3, ".gz") == 0) {
			file = gzopen(filename, "wx"); /* 'x' is for checking existance */
		} else {
			file = gzopen(filename, "wxT"); /* without compression */
		}
	}
	if (!file) {
		if (errno == EEXIST) {
			fprintf(stderr, "Error: Output file '%s' has already existed!\n", filename);
		} else {
			fprintf(stderr, "Error: Can not open output file '%s'\n", filename);
		}
		return 1;
	}

	if (output_cmap) {
		nick_map_write_cmap(file, map);
	} else {
		nick_map_write(file, map);
	}
	gzclose(file);
	return 0;
}
