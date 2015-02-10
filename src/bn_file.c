#include <ctype.h>
#include <assert.h>
#include "bn_file.h"
#include "version.h"

static inline int to_integer(double x) { return (int)(x + .5); }

static inline int string_begins_as(const char *s, const char *prefix)
{
	return (memcmp(s, prefix, strlen(prefix)) == 0);
}

static int bn_read_txt(struct file *fp, struct fragment *f)
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
		file_error(fp, "Unexpected EOF when read fragment number");
		return -EINVAL;
	}
	if (count > 0) {
		if (array_reserve(f->nicks, count)) {
			return -ENOMEM;
		}
		f->nicks.size = count - 1;
		for (i = 0, pos = 0; i < count; ++i) {
			if (read_double(fp, &value)) {
				file_error(fp, "Unexpected EOF when read size of fragment[%d]", i);
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

static int bn_read_tsv_header(struct file *fp, struct nick_map *map)
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

static int bn_read_tsv(struct file *fp, struct fragment *f)
{
	char name[sizeof(f->name)];
	char strandText[4];
	int c, label, pos, strand;

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
				file_error(fp, "Missing fragment end line");
				return -EINVAL;
			}
		}

		if (read_integer(fp, &label)) {
			file_error(fp, "Failed to read 'label' column");
			return -EINVAL;
		}

		if (read_integer(fp, &pos)) {
			file_error(fp, "Failed to read 'pos' column");
			return -EINVAL;
		}

		if (read_string(fp, strandText, sizeof(strandText))) {
			file_error(fp, "Failed to read 'strand' column");
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
			file_error(fp, "Unknown strand text '%s'", strandText);
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

static int bn_read_bnx(struct file *fp, struct fragment *f)
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
				file_error(fp, "Missing label info line");
				return -EINVAL;
			}
			if (read_string(fp, f->name, sizeof(f->name))) {
				file_error(fp, "Failed to read molecule ID");
				return -EINVAL;
			}
			if (read_double(fp, &value)) {
				file_error(fp, "Failed to read molecule size");
				return -EINVAL;
			}
			f->size = to_integer(value);
		} else if (strcmp(type, "1") == 0) { /* label positions */
			if (!f->name[0]) {
				file_error(fp, "Missing molecule info line");
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

static int bn_read_cmap_header(struct file *fp, struct nick_map *map)
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

static int bn_read_cmap(struct file *fp, struct fragment *f)
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
				file_error(fp, "Missing fragment end line");
				return -EINVAL;
			}
		}

		for (i = 0, channel = 0, pos = 0; i < 5; ++i) {
			if (read_integer(fp, &value)) {
				file_error(fp, "Failed to read data");
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

int bn_skip_comment_lines(struct file *fp)
{
	char buf[256];
	int c;
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c != '#') {
			gzungetc(c, fp->file);
			break;
		}
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	return 0;
}

int bn_read_header(struct file *fp, int *format, struct nick_map *map)
{
	char buf[256];

	*format = FORMAT_UNKNOWN;
	while (!gzeof(fp->file) && *format == FORMAT_UNKNOWN) {
		if (current_char(fp) != '#') {
			break;
		}
		if (read_line(fp, buf, sizeof(buf))) {
			break;
		}
		if (string_begins_as(buf, "##fileformat=MAPv0.1")) {
			*format = FORMAT_TSV;
		} else if (string_begins_as(buf, "# BNX File Version:")) {
			*format = FORMAT_BNX;
		} else if (string_begins_as(buf, "# CMAP File Version:")) {
			*format = FORMAT_CMAP;
		}
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	if (*format == FORMAT_UNKNOWN) {
		*format = FORMAT_TXT; /* try simple text format for unrecognized file */
	}

	switch (*format) {
	case FORMAT_TSV: return bn_read_tsv_header(fp, map);
	case FORMAT_CMAP: return bn_read_cmap_header(fp, map);
	case FORMAT_TXT:
	case FORMAT_BNX:
	default: return bn_skip_comment_lines(fp);
	}
}

int bn_read(struct file *fp, int format, struct fragment *f)
{
	assert(fp != NULL);
	assert(f != NULL);

	switch (format) {
	case FORMAT_TXT: return bn_read_txt(fp, f);
	case FORMAT_TSV: return bn_read_tsv(fp, f);
	case FORMAT_BNX: return bn_read_bnx(fp, f);
	case FORMAT_CMAP: return bn_read_cmap(fp, f);
	default: assert(0); return -1;
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

int nick_map_load(struct nick_map *map, const char *filename)
{
	struct file *fp;
	struct fragment fragment = { };
	int err, format;

	assert(map != NULL);
	assert(map->fragments.data == NULL);
	assert(map->fragments.size == 0);

	fp = file_open(filename);
	if (!fp) {
		return -1;
	}
	if ((err = bn_read_header(fp, &format, map)) != 0) {
		file_close(fp);
		return err;
	}
	while (bn_read(fp, format, &fragment) == 0) {
		if (array_reserve(map->fragments, map->fragments.size + 1)) {
			return -ENOMEM;
		}
		/* move data from 'fragment' to 'map->fragments.data' */
		memcpy(&map->fragments.data[map->fragments.size++],
				&fragment, sizeof(struct fragment));
		fragment.nicks.data = NULL;
		fragment.nicks.size = 0;
		fragment.nicks.capacity = 0;
	}
	file_close(fp);
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
	gzprintf(file, "#name\tlabel\tpos\tstrand\tsize\n");
	return 0;
}

static int save_fragment_as_tsv(gzFile file, const struct fragment *fragment)
{
	const char * const STRAND[] = { "?", "+", "-", "+/-" };
	size_t i;
	const struct nick *n;
	for (i = 0, n = NULL; i < fragment->nicks.size; ++i) {
		n = &fragment->nicks.data[i];
		gzprintf(file, "%s\t%zd\t%d\t%s\t%d\n",
				fragment->name, i, n->pos, STRAND[n->flag & 3],
				n->pos - (i == 0 ? 0 : (n - 1)->pos));
	}
	gzprintf(file, "%s\t%zd\t%d\t*\t%d\n",
			fragment->name, fragment->nicks.size, fragment->size, fragment->size - (n ? n->pos : 0));
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
