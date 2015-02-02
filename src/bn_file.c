#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include "bn_file.h"
#include "version.h"

static inline int to_integer(double x) { return (int)(x + .5); }

static inline int string_begins_as(const char *s, const char *prefix)
{
	return (memcmp(s, prefix, strlen(prefix)) == 0);
}

static void skip_spaces(struct bn_file *fp)
{
	int c;
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '\n') {
			++fp->line;
		} else if (!isspace(c)) {
			gzungetc(c, fp->file);
			break;
		}
	}
}

static void skip_current_line(struct bn_file *fp)
{
	int c;
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '\n') {
			++fp->line;
			break;
		}
	}
}

static int read_line(struct bn_file *fp, char *buf, size_t bufsize)
{
	size_t i = 0;
	int c;

	assert(fp != NULL);
	assert(buf != NULL);
	assert(bufsize > 1);

	while ((c = gzgetc(fp->file)) != EOF) {
		buf[i++] = (char)c;
		if (c == '\n') {
			++fp->line;
		}
		if (i + 1 == bufsize || c == '\n') {
			break;
		}
	}
	buf[i] = '\0';
	return (i > 0 ? 0 : -1);
}

static void skip_to_next_line(struct bn_file *fp, char *buf, size_t bufsize)
{
	while (strchr(buf, '\n') == NULL) {
		if (read_line(fp, buf, bufsize)) break;
	}
}

static int read_string(struct bn_file *fp, char *buf, size_t bufsize)
{
	size_t i = 0;
	int c;

	assert(fp != NULL);
	assert(buf != NULL);
	assert(bufsize > 0);

	skip_spaces(fp);
	while ((c = gzgetc(fp->file)) != EOF) {
		if (isspace(c)) {
			gzungetc(c, fp->file);
			break;
		} else if (i + 1 < bufsize) {
			buf[i++] = (char)c;
		}
	}
	buf[i] = '\0';
	return (i > 0 ? 0 : -1);
}

static int read_integer(struct bn_file *fp, int *value)
{
	int c;

	assert(fp != NULL);
	assert(value != NULL);

	skip_spaces(fp);
	c = gzgetc(fp->file);
	if (!isdigit(c)) {
		gzungetc(c, fp->file);
		return -1;
	}
	*value = c - '0';
	while ((c = gzgetc(fp->file)) != EOF) {
		if (!isdigit(c)) {
			gzungetc(c, fp->file);
			break;
		}
		*value = *value * 10 + (c - '0');
	}
	return 0;
}

static int read_double(struct bn_file *fp, double *value)
{
	int c, point = 0;
	double factor = 0.1;

	assert(fp != NULL);
	assert(value != NULL);

	skip_spaces(fp);
	c = gzgetc(fp->file);
	if (c != '.' && !isdigit(c)) {
		gzungetc(c, fp->file);
		return -1;
	}
	if (c == '.') {
		point = 1;
		*value = 0;
	} else {
		*value = c - '0';
	}
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '.') {
			if (point) {
				return -1;
			} else {
				point = 1;
			}
		} else if (isdigit(c)) {
			if (point) {
				*value += (c - '0') * factor;
				factor /= 10;
			} else {
				*value = *value * 10 + (c - '0');
			}
		} else if (isspace(c)) {
			gzungetc(c, fp->file);
			break;
		} else {
			return -1;
		}
	}
	return 0;
}

#define bn_file_error(fp, fmt, args...) \
	fprintf(stderr, "Error: " fmt " at line %zd of file '%s'\n", ##args, (fp)->line, (fp)->name)

struct bn_file *bn_open(const char *filename)
{
	struct bn_file *fp;
	int c;

	fp = malloc(sizeof(struct bn_file));
	if (!fp) {
		return NULL;
	}

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		fp->file = gzdopen(0, "r"); /* stdin */
	} else {
		fp->file = gzopen(filename, "r");
	}
	if (!fp->file) {
		fprintf(stderr, "Error: Can not open file '%s'!\n", filename);
		free(fp);
		return NULL;
	}

	fp->format = FORMAT_UNKNOWN;
	fp->name = filename;
	fp->line = 1;
	c = gzgetc(fp->file);
	gzungetc(c, fp->file);

	if (c != '#') { /* no any comment line, try as the simple text format */
		fp->format = FORMAT_TXT;
	} else {
		char buf[256];
		while (!gzeof(fp->file)) {
			if (read_line(fp, buf, sizeof(buf))) break;
			if (string_begins_as(buf, "##fileformat=MAPv0.1")) {
				fp->format = FORMAT_TSV;
			} else if (string_begins_as(buf, "# BNX File Version:")) {
				fp->format = FORMAT_BNX;
			} else if (string_begins_as(buf, "# CMAP File Version:")) {
				fp->format = FORMAT_CMAP;
			}
			skip_to_next_line(fp, buf, sizeof(buf));
			if (fp->format != FORMAT_UNKNOWN) break;
			c = gzgetc(fp->file);
			gzungetc(c, fp->file);
			if (c != '#') break;
		}
		if (fp->format == FORMAT_UNKNOWN) {
			fprintf(stderr, "Error: Unknown file format of '%s'!\n", filename);
			gzclose(fp->file);
			free(fp);
			return NULL;
		}
	}
	return fp;
}

static int bn_read_txt(struct bn_file *fp, struct fragment *f)
{
	/* format as: id, #intervals, size[1], ..., size[#intervals] */
	int count, i, pos;
	double value;

	f->name[0] = '\0';
	f->nicks.size = 0;
	f->size = 0;

	if (read_string(fp, f->name, sizeof(f->name))) {
		return -EINVAL;
	}
	if (read_integer(fp, &count)) {
		bn_file_error(fp, "Unexpected EOF when read fragment number");
		return -EINVAL;
	}
	if (count > 0) {
		if (array_reserve(f->nicks, count)) {
			return -ENOMEM;
		}
		f->nicks.size = count - 1;
		for (i = 0, pos = 0; i < count; ++i) {
			if (read_double(fp, &value)) {
				bn_file_error(fp, "Unexpected EOF when read size of fragment[%d]", i);
				return -EINVAL;
			}
			pos += to_integer(value);
			if (i + 1 < count) {
				f->nicks.data[i].pos = pos;
				f->nicks.data[i].flag = 0;
			} else {
				f->size = pos;
			}
		}
	}
	return 0;
}

static int bn_read_tsv_header(struct bn_file *fp, struct nick_map *map)
{
	char buf[256];
	int c;

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c != '#') {
			gzungetc(c, fp->file);
			break;
		}
		if (read_line(fp, buf, sizeof(buf))) break;
		if (string_begins_as(buf, "#enzyme=")) {
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
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	return 0;
}

static int bn_read_tsv(struct bn_file *fp, struct fragment *f)
{
	char name[sizeof(f->name)];
	char strandText[4];
	int c, pos, strand;

	f->name[0] = '\0';
	f->nicks.size = 0;
	f->size = 0;

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '#') {
			skip_current_line(fp);
			continue;
		}
		gzungetc(c, fp->file);

		if (read_string(fp, name, sizeof(name))) break;
		if (!f->name[0]) {
			snprintf(f->name, sizeof(f->name), "%s", name);
		} else {
			if (strncmp(f->name, name, sizeof(f->name) - 1) != 0) {
				bn_file_error(fp, "Missing fragment end line");
				return -EINVAL;
			}
		}

		if (read_integer(fp, &pos)) {
			bn_file_error(fp, "Failed to read pos column");
			return -EINVAL;
		}

		if (read_string(fp, strandText, sizeof(strandText))) {
			bn_file_error(fp, "Failed to read strand column");
		}
		if (strcmp(strandText, "?") == 0) {
			strand = 0;
		} else if (strcmp(strandText, "+") == 0) {
			strand = NICK_PLUS_STRAND;
		} else if (strcmp(strandText, "-") == 0) {
			strand = NICK_MINUS_STRAND;
		} else if (strcmp(strandText, "+/-") == 0) {
			strand = NICK_PLUS_STRAND | NICK_MINUS_STRAND;
		} else if (strcmp(strandText, "*") == 0) {
			strand = -1;
		} else {
			bn_file_error(fp, "Unknown strand text '%s'", strandText);
			return -EINVAL;
		}

		skip_current_line(fp);

		if (strand >= 0) {
			if (array_reserve(f->nicks, f->nicks.size + 1)) {
				return -ENOMEM;
			}
			f->nicks.data[f->nicks.size].pos = pos;
			f->nicks.data[f->nicks.size].flag = strand;
			++f->nicks.size;
		} else {
			f->size = pos;
			break;
		}
	}
	return (f->name[0] ? 0 : -1);
}

static int bn_read_bnx(struct bn_file *fp, struct fragment *f)
{
	char type[5];
	int c;
	double value;

	f->name[0] = '\0';
	f->nicks.size = 0;
	f->size = 0;

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '#') {
			skip_current_line(fp);
			continue;
		}
		gzungetc(c, fp->file);

		if (read_string(fp, type, sizeof(type))) break;
		if (strcmp(type, "0") == 0) { /* molecule info */
			if (f->name[0]) {
				bn_file_error(fp, "Missing label info line");
				return -EINVAL;
			}
			if (read_string(fp, f->name, sizeof(f->name))) {
				bn_file_error(fp, "Failed to read molecule ID");
				return -EINVAL;
			}
			if (read_double(fp, &value)) {
				bn_file_error(fp, "Failed to read molecule size");
				return -EINVAL;
			}
			f->size = to_integer(value);
		} else if (strcmp(type, "1") == 0) { /* label positions */
			if (!f->name[0]) {
				bn_file_error(fp, "Missing molecule info line");
				return -EINVAL;
			}
			while (read_double(fp, &value) == 0) {
				int pos = to_integer(value);
				if (pos == f->size) {
					break;
				}
				if (array_reserve(f->nicks, f->nicks.size + 1)) {
					return -ENOMEM;
				}
				f->nicks.data[f->nicks.size].pos = pos;
				f->nicks.data[f->nicks.size].flag = 0;
				++f->nicks.size;
			}
			skip_current_line(fp);
			break;
		}
		skip_current_line(fp);
	}
	return (f->name[0] ? 0 : -1);
}

static int bn_read_cmap_header(struct bn_file *fp, struct nick_map *map)
{
	char buf[256];
	int c;
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c != '#') {
			gzungetc(c, fp->file);
			break;
		}
		if (read_line(fp, buf, sizeof(buf))) break;
		if (string_begins_as(buf, " Nickase Recognition Site 1:")) {
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
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	return 0;
}

static int bn_read_cmap(struct bn_file *fp, struct fragment *f)
{
	char map_id[sizeof(f->name)];
	int c, i, value, channel, pos;

	f->name[0] = '\0';
	f->nicks.size = 0;
	f->size = 0;

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '#') {
			skip_current_line(fp);
			continue;
		}
		gzungetc(c, fp->file);

		if (read_string(fp, map_id, sizeof(map_id))) break;
		if (!f->name[0]) {
			snprintf(f->name, sizeof(f->name), "%s", map_id);
		} else {
			if (strncmp(f->name, map_id, sizeof(f->name) - 1) != 0) {
				bn_file_error(fp, "Missing fragment end line");
				return -EINVAL;
			}
		}

		for (i = 0, channel = 0, pos = 0; i < 5; ++i) {
			if (read_integer(fp, &value)) {
				bn_file_error(fp, "Failed to read data");
				return -EINVAL;
			}
			switch (i) {
			case 3: channel = value; break;
			case 4: pos = value; break;
			default: break;
			}
		}
		if (channel == 1) {
			if (array_reserve(f->nicks, f->nicks.size + 1)) {
				return -ENOMEM;
			}
			f->nicks.data[f->nicks.size].pos = pos;
			f->nicks.data[f->nicks.size].flag = 0;
			++f->nicks.size;
		} else {
			assert(channel == 0);
			f->size = pos;
			skip_current_line(fp);
			break;
		}
		skip_current_line(fp);
	}
	return (f->name[0] ? 0 : -1);
}

int bn_read_header(struct bn_file *fp, struct nick_map *map)
{
	switch (fp->format) {
	case FORMAT_TXT: return 0;
	case FORMAT_TSV: return bn_read_tsv_header(fp, map);
	case FORMAT_BNX: return 0;
	case FORMAT_CMAP: return bn_read_cmap_header(fp, map);
	default: assert(0); return -1;
	}
}

int bn_read(struct bn_file *fp, struct fragment *f)
{
	assert(fp != NULL);
	assert(f != NULL);

	switch (fp->format) {
	case FORMAT_TXT: return bn_read_txt(fp, f);
	case FORMAT_TSV: return bn_read_tsv(fp, f);
	case FORMAT_BNX: return bn_read_bnx(fp, f);
	case FORMAT_CMAP: return bn_read_cmap(fp, f);
	default: assert(0); return -1;
	}
}

void bn_close(struct bn_file *fp)
{
	if (fp) {
		if (fp->file) {
			gzclose(fp->file);
		}
		free(fp);
	}
}

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

int load_name_list(struct name_list *name_list, const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		return -EINVAL;
	}
	while (!feof(fp)) {
		char buf[256];
		size_t len;
		if (!fgets(buf, sizeof(buf), fp)) break;
		if (buf[0] == '#') continue;
		len = strlen(buf);
		if (len > 0 && buf[len - 1] == '\n') {
			buf[len - 1] = '\0';
		}
		if (array_reserve(name_list->names, name_list->names.size + 1)) {
			fclose(fp);
			return -ENOMEM;
		}
		name_list->names.data[name_list->names.size] = strdup(buf);
		++name_list->names.size;
	}
	fclose(fp);
	return 0;
}

void free_name_list(struct name_list *name_list)
{
	size_t i;
	for (i = 0; i < name_list->names.size; ++i) {
		free(name_list->names.data[i]);
	}
	array_free(name_list->names);
}

int name_list_has(const struct name_list *name_list, const char *name)
{
	size_t i;
	for (i = 0; i < name_list->names.size; ++i) {
		if (strcmp(name, name_list->names.data[i]) == 0) {
			return 1;
		}
	}
	return 0;
}

int nick_map_load(struct nick_map *map, const char *filename)
{
	struct bn_file *fp;
	struct fragment fragment = { };
	int err;

	fp = bn_open(filename);
	if (!fp) {
		return -1;
	}
	if ((err = bn_read_header(fp, map)) != 0) {
		bn_close(fp);
		return err;
	}
	while (bn_read(fp, &fragment) == 0) {
		if (array_reserve(map->fragments, map->fragments.size + 1)) {
			return -ENOMEM;
		}
		memcpy(&map->fragments.data[map->fragments.size++], &fragment, sizeof(struct fragment));
	}
	bn_close(fp);
	return 0;
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

static int save_fragment_as_txt(gzFile file, const struct fragment *fragment)
{
	const struct nick *n;
	size_t i;
	gzprintf(file, "%s %d", fragment->name, fragment->nicks.size + 1);
	for (i = 0, n = NULL; i < fragment->nicks.size; ++i) {
		n = &fragment->nicks.data[i];
		gzprintf(file, " %d", n->pos - (i == 0 ? 0 : (n - 1)->pos));
	}
	gzprintf(file, " %d\n", fragment->size - (n ? n->pos : 0));
	return 0;
}

static int save_as_txt(gzFile file, const struct nick_map *map)
{
	size_t i;
	for (i = 0; i < map->fragments.size; ++i) {
		save_fragment_as_txt(file, &map->fragments.data[i]);
	}
	return 0;
}

static int save_tsv_header(gzFile file, const struct nick_map *map)
{
	gzprintf(file, "##fileformat=MAPv0.1\n");
	if (map->enzyme[0] && map->rec_seq[0]) {
		gzprintf(file, "##enzyme=%s/%s\n", map->enzyme, map->rec_seq);
	}
	gzprintf(file, "##program=bntools\n");
	gzprintf(file, "##programversion="VERSION"\n");
	write_command_line(file);
	gzprintf(file, "#name\tpos\tstrand\tsize\n");
	return 0;
}

static int save_fragment_as_tsv(gzFile file, const struct fragment *fragment)
{
	const char * const STRAND[] = { "?", "+", "-", "+/-" };
	size_t i;
	const struct nick *n;
	for (i = 0, n = NULL; i < fragment->nicks.size; ++i) {
		n = &fragment->nicks.data[i];
		gzprintf(file, "%s\t%d\t%s\t%d\n",
				fragment->name, n->pos, STRAND[n->flag & 3],
				n->pos - (i == 0 ? 0 : (n - 1)->pos));
	}
	gzprintf(file, "%s\t%d\t*\t%d\n",
			fragment->name, fragment->size, fragment->size - (n ? n->pos : 0));
	return 0;
}

static int save_as_tsv(gzFile file, const struct nick_map *map)
{
	size_t i;
	save_tsv_header(file, map);
	for (i = 0; i < map->fragments.size; ++i) {
		save_fragment_as_tsv(file, &map->fragments.data[i]);
	}
	return 0;
}

static int save_bnx_header(gzFile file, const struct nick_map *map)
{
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
	return 0;
}

static int save_fragment_as_bnx(gzFile file, const struct fragment *fragment)
{
	size_t i;
	gzprintf(file, "0\t%s\t%zd\n1", fragment->name, fragment->size);
	for (i = 0; i < fragment->nicks.size; ++i) {
		gzprintf(file, "\t%d", fragment->nicks.data[i].pos);
	}
	gzprintf(file, "\t%d\n", fragment->size);
	return 0;
}

static int save_as_bnx(gzFile file, const struct nick_map *map)
{
	size_t i;
	save_bnx_header(file, map);
	for (i = 0; i < map->fragments.size; ++i) {
		save_fragment_as_bnx(file, &map->fragments.data[i]);
	}
	return 0;
}

static int save_cmap_header(gzFile file, const struct nick_map *map)
{
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
	return 0;
}

static int save_fragment_as_cmap(gzFile file, const struct fragment *fragment)
{
	size_t i;
	for (i = 0; i < fragment->nicks.size; ++i) {
		gzprintf(file, "%s\t%d\t%zd\t%zd\t%d\t%d\t%d\t%d\t%d\n",
				fragment->name, fragment->size, fragment->nicks.size,
				i + 1, 1, fragment->nicks.data[i].pos, 0, 0, 0);
	}
	gzprintf(file, "%s\t%d\t%zd\t%zd\t%d\t%d\t%d\t%d\t%d\n",
			fragment->name, fragment->size, fragment->nicks.size,
			fragment->nicks.size + 1, 0, fragment->size, 0, 1, 1);
	return 0;
}

static int save_as_cmap(gzFile file, const struct nick_map *map)
{
	size_t i;
	save_cmap_header(file, map);
	for (i = 0; i < map->fragments.size; ++i) {
		save_fragment_as_cmap(file, &map->fragments.data[i]);
	}
	return 0;
}

gzFile open_gzfile_write(const char *filename)
{
	gzFile file;

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdout") == 0) {
		file = gzdopen(1, "wT"); /* 1: stdout; 'T': without compression */
	} else {
		size_t len = strlen(filename);
		if (len > 3 && strcmp(filename + len - 3, ".gz") == 0) {
			file = gzopen(filename, "wx9"); /* 'x': check existed; '9': best compression */
		} else {
			file = gzopen(filename, "wxT"); /* 'T': without compression */
		}
	}
	if (!file) {
		if (errno == EEXIST) {
			fprintf(stderr, "Error: Output file '%s' has already existed!\n", filename);
		} else {
			fprintf(stderr, "Error: Can not open output file '%s'\n", filename);
		}
		return NULL;
	}
	return file;
}

int nick_map_save(const struct nick_map *map, const char *filename, int format)
{
	gzFile file;
	int ret;

	file = open_gzfile_write(filename);
	if (!file) {
		return -1;
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

int save_header(gzFile file, const struct nick_map *map, int format)
{
	switch (format) {
	case FORMAT_TXT: return 0;
	case FORMAT_TSV: return save_tsv_header(file, map);
	case FORMAT_BNX: return save_bnx_header(file, map);
	case FORMAT_CMAP: return save_cmap_header(file, map);
	default: assert(0); return -EINVAL;
	}
}

int save_fragment(gzFile file, const struct fragment *fragment, int format)
{
	switch (format) {
	case FORMAT_TXT: return save_fragment_as_txt(file, fragment);
	case FORMAT_TSV: return save_fragment_as_tsv(file, fragment);
	case FORMAT_BNX: return save_fragment_as_bnx(file, fragment);
	case FORMAT_CMAP: return save_fragment_as_cmap(file, fragment);
	default: assert(0); return -EINVAL;
	}
}
