#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include <zlib.h>
#include "nick_map.h"
#include "version.h"

int parse_format_text(const char *s)
{
	if (strcmp(s, "txt") == 0) {
		return FORMAT_TXT;
	} else if (strcmp(s, "tsv") == 0) {
		return FORMAT_TSV;
	} else if (strcmp(s, "bnx") == 0) {
		return FORMAT_BNX;
	} else if (strcmp(s, "cmap") == 0) {
		return FORMAT_CMAP;
	} else {
		return FORMAT_UNKNOWN;
	}
}

static inline int to_integer(double x) { return (int)(x + .5); }

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

static int load_txt(gzFile file, struct nick_map *map)
{
	/* format as: id, #intervals, size[1], ..., size[#intervals] */
	struct fragment *f = NULL;
	char name[32] = "";
	double sum = 0;
	double value;
	int count = 0, i = 0, err, c;
	char buf[256];
	size_t pos = 0;

	while ((c = gzgetc(file)) != EOF) {
		if (!isspace(c)) {
			if (pos >= sizeof(buf)) {
				fprintf(stderr, "Error: Unexpected long word!\n");
				return -EINVAL;
			}
			buf[pos++] = (char)c;
		} else if (pos > 0) {
			buf[pos] = '\0';
			pos = 0;
			if (!name[0]) {
				snprintf(name, sizeof(name), "%s", buf);
				count = 0;
				i = 0;
				f = nick_map_add_fragment(map, name);
				if (!f) {
					return -ENOMEM;
				}
			} else if (count == 0) {
				count = atoi(buf);
				i = 0;
				if (count <= 0) {
					name[0] = '\0';
					count = 0;
				}
			} else {
				assert(f != NULL);
				value = atof(buf);
				sum += value;
				++i;
				if (i < count) {
					if ((err = nick_map_add_site(f,
							to_integer(sum), STRAND_UNKNOWN)) != 0) {
						return err;
					}
				} else {
					if ((err = nick_map_add_site(f,
							to_integer(sum), STRAND_END)) != 0) {
						return err;
					}
					f->size = to_integer(sum);
					sum = 0;
					name[0] = '\0';
					count = 0;
					i = 0;
				}
			}
		}
	}
	if (name[0] || i < count) {
		fprintf(stderr, "Error: Unexpected EOF!\n");
		return -EINVAL;
	}
	return 0;
}

static int load_tsv(gzFile file, long long lineNo, struct nick_map *map)
{
	char buf[256];
	struct fragment *f = NULL;
	char lastName[32] = "";

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
					nick_map_set_enzyme(map, p, q);
				}
			}
			continue;
		} else {
			char name[32];
			int pos;
			char strandText[8];
			int strand;

			if (sscanf(buf, "%s%d%s", name, &pos, strandText) != 3) {
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

			if (strcmp(lastName, name) != 0) {
				f = nick_map_add_fragment(map, name);
				snprintf(lastName, sizeof(lastName), "%s", name);
			}
			assert(f != NULL);

			nick_map_add_site(f, pos, strand);
			if (strand == STRAND_END) {
				f->size = pos;
			}
		}
	}
	return 0;
}

static int load_bnx(gzFile file, long long lineNo, struct nick_map *map)
{
	char buf[256];
	char molecule_id[64] = "";
	double length = 0;
	int molecule_length;
	struct fragment *f = NULL;
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
			f = nick_map_add_fragment(map, molecule_id);
			f->size = molecule_length;
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
					p = buf;
					continue;
				}
				if (sscanf(p, "%lf", &value) != 1) {
					fprintf(stderr, "Error: Failed in reading float value on line %lld\n", lineNo);
					return 1;
				}
				nick_map_add_site(f, to_integer(value),
						(*q == '\t' ? STRAND_UNKNOWN : STRAND_END));
				if (*q == '\t') {
					p = q + 1;
				} else {
					assert(*q == '\n');
					f->size = to_integer(value);
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
	struct fragment *f = NULL;
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
					nick_map_set_enzyme(map, p, q);
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
				f = nick_map_add_fragment(map, mapId);
			}
			assert(f != NULL);

			assert(labelChannel == 0 || labelChannel == 1);
			nick_map_add_site(f, position, (labelChannel == 1 ? 0 : STRAND_END));
			if (labelChannel == 0) {
				f->size = position;
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
	int c;

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		file = gzdopen(0, "r"); /* stdin */
	} else {
		file = gzopen(filename, "r");
	}
	if (!file) {
		fprintf(stderr, "Error: Can not open file '%s'\n", filename);
		return 1;
	}

	c = gzgetc(file);
	gzungetc(c, file);
	if (c != '#') {  /* no any comment line, try as the simple format */
		ret = load_txt(file, map);
	} else {
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
	}
	if (ret) {
		nick_map_free(map);
	}
	gzclose(file);
	return ret;
}

static void write_command_line(gzFile file)
{
	char name[64];
	FILE *fp;
	snprintf(name, sizeof(name), "/proc/%d/cmdline", getpid());
	fp = fopen(name, "r");
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

static int save_as_txt(gzFile file, const struct nick_map *map)
{
	size_t i, j;
	for (i = 0; i < map->fragments.size; ++i) {
		const struct fragment *f = &map->fragments.data[i];
		assert(f->nicks.size > 1);
		assert(f->nicks.data[0].pos == 0);
		assert(f->nicks.data[f->nicks.size - 1].pos == f->size);
		gzprintf(file, "%s %d", f->name, f->nicks.size - 1);
		for (j = 1; j < f->nicks.size; ++j) {
			const struct nick *n = f->nicks.data;
			gzprintf(file, " %d", n[j].pos - (j == 0 ? 0 : n[j - 1].pos));
		}
		gzprintf(file, "\n");
	}
	return 0;
}

static int save_as_tsv(gzFile file, const struct nick_map *map)
{
	static const char * const STRAND[] = { "?", "+", "-", "+/-", "*" };
	size_t i, j;

	gzprintf(file, "##fileformat=MAPv0.1\n");
	if (map->enzyme[0] && map->rec_seq[0]) {
		gzprintf(file, "##enzyme=%s/%s\n", map->enzyme, map->rec_seq);
	}
	gzprintf(file, "##program=bntools\n");
	gzprintf(file, "##programversion="VERSION"\n");
	write_command_line(file);
	gzprintf(file, "#name\tpos\tstrand\tsize\n");

	for (i = 0; i < map->fragments.size; ++i) {
		const struct fragment *f = &map->fragments.data[i];
		assert(f->nicks.size > 1);
		assert(f->nicks.data[0].pos == 0);
		assert(f->nicks.data[f->nicks.size - 1].pos == f->size);
		for (j = 1; j < f->nicks.size; ++j) {
			const struct nick *n = f->nicks.data;
			gzprintf(file, "%s\t%d\t%s\t%d\n",
					f->name, n[j].pos, STRAND[n[j].strand], n[j].pos - n[j - 1].pos);
		}
	}
	return 0;
}

static int save_as_bnx(gzFile file, const struct nick_map *map)
{
	size_t i, j;

	gzprintf(file, "# BNX File Version: 0.1\n");
	gzprintf(file, "# Label Channels: 1\n");
	if (map->enzyme[0] && map->rec_seq[0]) {
		gzprintf(file, "# Nickase Recognition Site 1: %s/%s\n", map->enzyme, map->rec_seq);
	} else {
		gzprintf(file, "# Nickase Recognition Site 1: unknown\n");
	}
	gzprintf(file, "# Number of Nanomaps: %zd\n", map->fragments.size);
	gzprintf(file, "#0h\tLabel Channel\tMapID\tLength\n");
	gzprintf(file, "#0f\tint\tint\tfloat\n");
	gzprintf(file, "#1h\tLabel Channel\tLabelPositions[N]\n");
	gzprintf(file, "#1f\tint\tfloat\n");

	for (i = 0; i < map->fragments.size; ++i) {
		const struct fragment *f = &map->fragments.data[i];
		assert(f->nicks.size > 1);
		assert(f->nicks.data[0].pos == 0);
		assert(f->nicks.data[f->nicks.size - 1].pos == f->size);
		gzprintf(file, "0\t%s\t%zd\n1", f->name, f->size - 1);
		for (j = 1; j < f->nicks.size; ++j) {
			gzprintf(file, "\t%d", f->nicks.data[j].pos);
		}
		gzprintf(file, "\n");
	}
	return 0;
}

static int save_as_cmap(gzFile file, const struct nick_map *map)
{
	size_t i, j;

	gzprintf(file, "# CMAP File Version:  0.1\n");
	gzprintf(file, "# Label Channels:  1\n");
	if (map->enzyme[0] && map->rec_seq[0]) {
		gzprintf(file, "# Nickase Recognition Site 1:  %s/%s\n", map->enzyme, map->rec_seq);
	} else {
		gzprintf(file, "# Nickase Recognition Site 1:  unknown\n");
	}
	gzprintf(file, "# Number of Consensus Nanomaps:    %zd\n", map->fragments.size);
	gzprintf(file, "#h CMapId\tContigLength\tNumSites\tSiteID"
			"\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n");
	gzprintf(file, "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n");

	for (i = 0; i < map->fragments.size; ++i) {
		const struct fragment *f = &map->fragments.data[i];
		assert(f->nicks.size > 1);
		assert(f->nicks.data[0].pos == 0);
		assert(f->nicks.data[f->nicks.size - 1].pos == f->size);
		for (j = 1; j < f->nicks.size; ++j) {
			gzprintf(file, "%s\t%d\t%zd\t%zd\t%d\t%d\t%d\t%d\t%d\n",
					f->name,                          /* cmapId */
					f->size,                          /* contigLength */
					f->nicks.size - 1,                /* numSites */
					j,                                /* siteId */
					(j < f->nicks.size - 1 ? 1 : 0),  /* channel */
					f->nicks.data[j].pos,             /* position */
					0,                                /* stddev */
					(j < f->nicks.size - 1 ? 0 : 1),  /* coverage */
					(j < f->nicks.size - 1 ? 0 : 1)); /* occurance */
		}
	}
	return 0;
}

int nick_map_save(const struct nick_map *map, const char *filename, int format)
{
	gzFile file;
	int ret;

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

	switch (format) {
	case FORMAT_TXT: ret = save_as_txt(file, map); break;
	case FORMAT_TSV: ret = save_as_tsv(file, map); break;
	case FORMAT_BNX: ret = save_as_bnx(file, map); break;
	case FORMAT_CMAP: ret = save_as_cmap(file, map); break;
	default: assert(0); ret = -EINVAL; break;
	}
	gzclose(file);

	if (ret) {
		unlink(filename);
	}
	return ret;
}
