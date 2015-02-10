#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "version.h"
#include "ref_map.h"
#include "base_map.h"
#include "io_base.h"
#include "bn_file.h"

void ref_map_init(struct ref_map *ref)
{
	memset(ref, 0, sizeof(struct ref_map));
}

void ref_map_free(struct ref_map *ref)
{
	array_free(ref->index_);
	array_free(ref->nodes);
	nick_map_free(&ref->map);
}

int prepare_rec_site(struct rec_site *site, const char *enzyme, const char *rec_seq)
{
	int i;

	assert(site != NULL);
	assert(enzyme != NULL);
	assert(rec_seq != NULL);

	memset(site, 0, sizeof(struct rec_site));
	snprintf(site->enzyme, sizeof(site->enzyme), "%s", enzyme);
	snprintf(site->rec_seq, sizeof(site->rec_seq), "%s", rec_seq);
	assert(strcmp(site->enzyme, enzyme) == 0); /* ensure no buffer overflow */
	assert(strcmp(site->rec_seq, rec_seq) == 0);

	site->nick_offset = -1;
	site->rec_seq_size = 0;
	for (i = 0; rec_seq[i]; ++i) {
		if (rec_seq[i] == '^') {
			if (site->nick_offset >= 0) {
				fprintf(stderr, "Error: Invalid recognition sequence '%s'\n", rec_seq);
				return 1;
			}
			site->nick_offset = i;
		} else {
			char c = char_to_base(rec_seq[i]);
			if (c == 0) {
				fprintf(stderr, "Error: Invalid character '%c' in recognition sequence '%s'\n", rec_seq[i], rec_seq);
				return 1;
			}
			if (site->rec_seq_size >= sizeof(site->rec_bases)) {
				fprintf(stderr, "Error: Recognition sequence is too long\n");
				return 1;
			}
			site->rec_bases[site->rec_seq_size++] = c;
		}
	}
	if (site->nick_offset < 0) {
		fprintf(stderr, "Error: Missing '^' in recognition sequence '%s'\n", rec_seq);
		return 1;
	}
	for (i = 0, site->palindrome = 1; i < site->rec_seq_size / 2; ++i) {
		if (site->rec_bases[i] != base_to_comp(site->rec_bases[site->rec_seq_size - i - 1])) {
			site->palindrome = 0;
			break;
		}
	}
	return 0;
}

static int seq_match(const char *ref, const char *query, size_t len, int strand)
{
	size_t i;
	for (i = 0; i < len; ++i) {
		char r = ref[i];
		char q = (strand == 0 ? query[i] : base_to_comp(query[len - i - 1]));
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

static int is_chrom(const char *name)
{
	const char *p = name;
	if (memcmp(p, "chr", 3) == 0) {
		p += 3;
	}
	if (isdigit(p[0]) && p[1] == '\0') {
		return 1;
	} else if (isdigit(p[0]) && isdigit(p[1]) && p[2] == '\0') {
		return 1;
	} else if ((p[0] == 'X' || p[1] == 'Y') && p[2] == '\0') {
		return 1;
	} else {
		return 0;
	}
}

int nick_map_load_seq(struct ref_map *ref, const char *filename,
		const struct rec_site *site, int chrom_only, int verbose)
{
	struct file *fp;
	struct fragment *f = NULL;
	char name[MAX_CHROM_NAME_SIZE] = "";
	struct buffer buf = { };
	int c, ret = 0, base_count = 0;
	int format = 0; /* 1 - FASTA, 2 - FASTQ */

	if (!ref->map.enzyme[0]) {
		assert(sizeof(ref->map.enzyme) == MAX_ENZYME_NAME_SIZE + 1);
		assert(sizeof(ref->map.rec_seq) == MAX_REC_SEQ_SIZE + 1);
		assert(sizeof(site->enzyme) == MAX_ENZYME_NAME_SIZE + 1);
		assert(sizeof(site->rec_seq) == MAX_REC_SEQ_SIZE + 1);
		memcpy(ref->map.enzyme, site->enzyme, sizeof(ref->map.enzyme));
		memcpy(ref->map.rec_seq, site->rec_seq, sizeof(ref->map.rec_seq));
	}

	fp = file_open(filename);
	if (!fp) {
		return 1;
	}

	c = gzgetc(fp->file);
	if (c == '>') {
		format = 1;
	} else if (c == '@') {
		format = 2;
	} else {
		fprintf(stderr, "Error: File '%s' is not in FASTA/FASTQ format\n", filename);
		ret = -EINVAL;
		goto out;
	}

	for (;;) {
		int newline = 1;

		if (read_string(fp, name, sizeof(name))) {
			break;
		}
		skip_current_line(fp);

		if (chrom_only && !is_chrom(name)) {
			f = NULL;
		} else {
			f = nick_map_add_fragment(&ref->map, name);
			if (!f) {
				ret = -ENOMEM;
				goto out;
			}
			base_count = 0;
			if (verbose > 0) {
				fprintf(stderr, "Loading fragment '%s' ... ", name);
			}
		}

		for (;;) {
			int base;

			c = gzgetc(fp->file);
			if (c == EOF) {
				break;
			} else if (c == '\n') {
				++fp->line;
			}

			base = char_to_base(c);
			if (base) {
				if (buf.pos >= sizeof(buf.data)) {
					memcpy(buf.data, buf.data + sizeof(buf.data)
							- site->rec_seq_size + 1,
							site->rec_seq_size - 1);
					buf.pos = site->rec_seq_size - 1;
				}
				++base_count;
				buf.data[buf.pos++] = base;

				if (buf.pos >= site->rec_seq_size) {
					int strand, matched;
					for (strand = 0, matched = 0; strand < 2; ++strand) {
						if (matched || seq_match(buf.data + buf.pos - site->rec_seq_size,
								site->rec_bases, site->rec_seq_size, strand)) {
							int site_pos = base_count - (strand == 1 ? site->nick_offset
									: (site->rec_seq_size - site->nick_offset));
							if (nick_map_add_site(f, site_pos,
									(strand == 0 ? NICK_PLUS_STRAND : NICK_MINUS_STRAND))) {
								return -ENOMEM;
							}
							matched = site->palindrome;
						}
					}
				}
			} else {
				if (format == 1) {
					if (newline && c == '>') {
						break;
					}
				} else {
					assert(format == 2);
					if (c == '\n') {
						skip_current_line(fp); /* skip 3rd line */
						skip_current_line(fp); /* skip 4rd line */
						c = gzgetc(fp->file);
						break;
					}
				}
			}

			newline = (c == '\n');
		}

		f->size = base_count;
		if (f && verbose > 0) {
			fprintf(stderr, "%d bp\n", base_count);
		}
	}
out:
	file_close(fp);
	return ret;
}

static int meet_last(const struct ref_index *p, int i)
{
	if (p->direct > 0) {
		return (p->node[i].flag & LAST_INTERVAL) != 0;
	} else {
		assert(p->direct < 0);
		return (p->node[i].flag & FIRST_INTERVAL) != 0;
	}
}

static int sort_by_size(const void *a, const void *b)
{
	const struct ref_index *pa = a;
	const struct ref_index *pb = b;
	int i, j;
	for (i = 0, j = 0; ; i += pa->direct, j += pb->direct) {
		if (pa->node[i].size < pb->node[j].size) return -1;
		if (pa->node[i].size > pb->node[j].size) return 1;
		if (meet_last(pa, i) && meet_last(pb, j)) return 0;
		if (meet_last(pa, i)) return -1;
		if (meet_last(pb, j)) return 1;
	}
	return 0;
}

static void set_node(struct ref_node *node,
		size_t chrom, size_t label, int pos, int size, int flag)
{
	node->chrom = chrom;
	node->label = label;
	node->pos = pos;
	node->size = size;
	node->flag = flag;
}

static int ref_map_prepare_nodes(struct ref_map *ref)
{
	size_t count, i, j;

	assert(ref->nodes.size == 0);

	for (count = 0, i = 0; i < ref->map.fragments.size; ++i) {
		if (ref->map.fragments.data[i].nicks.size > 1) {
			count += ref->map.fragments.data[i].nicks.size + 1;
		}
	}
	if (array_reserve(ref->nodes, count)) {
		return -ENOMEM;
	}

	for (i = 0; i < ref->map.fragments.size; ++i) {
		const struct fragment *f = &ref->map.fragments.data[i];
		if (f->nicks.size <= 1) continue;
		set_node(&ref->nodes.data[ref->nodes.size++],
				i, 0, 0, f->nicks.data[0].pos, FIRST_INTERVAL);
		for (j = 0; j + 1 < f->nicks.size; ++j) {
			set_node(&ref->nodes.data[ref->nodes.size++],
					i, j + 1, f->nicks.data[j].pos,
					f->nicks.data[j + 1].pos - f->nicks.data[j].pos, 0);
		}
		set_node(&ref->nodes.data[ref->nodes.size++],
				i, f->nicks.size, f->nicks.data[f->nicks.size - 1].pos,
				f->size - f->nicks.data[f->nicks.size - 1].pos, LAST_INTERVAL);
	}
	assert(ref->nodes.size == count);
	return 0;
}

int ref_map_build_index(struct ref_map *ref)
{
	size_t count, i, j, k, m, n;

	if (ref_map_prepare_nodes(ref)) {
		return -ENOMEM;
	}

	assert(ref->index_.size == 0);

	count = 0;
	for (i = 0; i < ref->map.fragments.size; ++i) {
		if (ref->map.fragments.data[i].nicks.size > 1) {
			count += (ref->map.fragments.data[i].nicks.size - 1) * 2;
		}
	}
	if (array_reserve(ref->index_, count)) {
		return -ENOMEM;
	}

	for (i = 0, m = 0, n = 0; i < ref->map.fragments.size; ++i) {
		const struct fragment *f = &ref->map.fragments.data[i];
		if (f->nicks.size <= 1) continue;
		++m;
		for (j = 0; j + 1 < f->nicks.size; ++j) {
			for (k = 0; k < 2; ++k) {
				ref->index_.data[n].node = &ref->nodes.data[m];
				ref->index_.data[n].direct = (k == 0 ? 1 : -1);
				ref->index_.data[n].uniq_count = 0;
				++n;
			}
			++m;
		}
		++m;
	}
	assert(n == count);
	ref->index_.size = count;

	qsort(ref->index_.data, ref->index_.size, sizeof(struct ref_index), sort_by_size);

	for (i = 0; i + 1 < ref->index_.size; ++i) {
		struct ref_index *a = &ref->index_.data[i];
		struct ref_index *b = &ref->index_.data[i + 1];
		int x, y, z;
		for (x = 0, y = 0, z = 0; ;
				x += ref->index_.data[i].direct,
				y += ref->index_.data[i + 1].direct, ++z) {
			if (a->node[x].size != b->node[y].size) break;
			if (meet_last(a, x) || meet_last(b, y)) {
				++z;
				break;
			}
		}
		if (a->uniq_count < z + 1) {
			a->uniq_count = z + 1;
		}
		if (b->uniq_count < z + 1) {
			b->uniq_count = z + 1;
		}
	}
	return 0;
}

const char *get_index_filename(const char *filename, char *buf, size_t bufsize)
{
	struct stat sb;
	size_t len;

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		snprintf(buf, bufsize, "-");
	} else if (stat(filename, &sb) == 0 && !S_ISREG(sb.st_mode)) {
		snprintf(buf, bufsize, "-");
	} else {
		snprintf(buf, bufsize, "%s", filename);
		len = strlen(filename);
		if (len > 3 && strcmp(filename + len - 3, ".gz") == 0) {
			len -= 3;
		}
		snprintf(buf + len, bufsize - len, ".idx.gz");
	}
	return buf;
}

int ref_map_save(const struct ref_map *ref, const char *filename)
{
	gzFile file;
	size_t i, j;

	file = open_gzfile_write(filename);
	if (!file) {
		return -EINVAL;
	}

	gzprintf(file, "##fileformat=IDXv0.1\n");
	gzprintf(file, "##program=bntools\n");
	gzprintf(file, "##programversion="VERSION"\n");
	gzprintf(file, "#index\tchrom\tlabel\tstrand\tname\tpos\tsize\tuniq\tseq\n");

	for (i = 0; i < ref->index_.size; ++i) {
		const struct ref_index *r = &ref->index_.data[i];
		gzprintf(file, "%zd\t%zd\t%zd\t%s\t%s\t%d\t%d\t%d\t",
				r->node - ref->nodes.data,
				r->node->chrom + 1, r->node->label, (r->direct > 0 ? "+" : "-"),
				ref->map.fragments.data[r->node->chrom].name,
				r->node->pos, r->node->size, r->uniq_count);
		for (j = 0; j < r->uniq_count; ++j) {
			const struct ref_node *n = &r->node[j * r->direct];
			gzprintf(file, "%s%d", (j == 0 ? "": ","), n->size);
			if (meet_last(r, j * r->direct)) break;
		}
		gzprintf(file, "\n");
	}
	gzclose(file);
	return 0;
}

int ref_map_load(struct ref_map *ref, const char *filename)
{
	struct file *file;
	size_t count;
	int value;
	size_t index, chrom, label;
	char directText[2];
	int direct;
	char name[64];
	int pos, size, uniq;
	size_t i, m;
	const struct ref_node *node;

	if (ref_map_prepare_nodes(ref)) {
		return -ENOMEM;
	}

	assert(ref->index_.size == 0);

	count = 0;
	for (i = 0; i < ref->map.fragments.size; ++i) {
		if (ref->map.fragments.data[i].nicks.size > 1) {
			count += (ref->map.fragments.data[i].nicks.size - 1) * 2;
		}
	}
	if (array_reserve(ref->index_, count)) {
		return -ENOMEM;
	}

	file = file_open(filename);
	if (!file) {
		return -EINVAL;
	}
	if (bn_skip_comment_lines(file)) {
		file_close(file);
		return -EINVAL;
	}
	for (m = 0;;) {
		if (read_integer(file, &value)) {
			break;
		}
		if (value < 0) {
			file_error(file, "Invalid value in 'index' column");
			file_close(file);
			return -EINVAL;
		}
		index = value;
		node = ref->nodes.data + index;

		if (read_integer(file, &value)) {
			file_error(file, "Failed to read 'chrom' column");
			file_close(file);
			return -EINVAL;
		}
		if (value <= 0) {
			file_error(file, "Invalid value in 'chrom' column");
			file_close(file);
			return -EINVAL;
		}
		chrom = value - 1;
		if (node->chrom != chrom) {
			file_error(file, "Column 'chrom' does not match");
			file_close(file);
			return -EINVAL;
		}

		if (read_integer(file, &value)) {
			file_error(file, "Failed to read 'label' column");
			file_close(file);
			return -EINVAL;
		}
		label = value;
		if (node->label != label) {
			file_error(file, "Column 'label' does not match");
			file_close(file);
			return -EINVAL;
		}

		if (read_string(file, directText, sizeof(directText))) {
			file_error(file, "Failed to read 'strand' column");
			file_close(file);
			return -EINVAL;
		}
		if (strcmp(directText, "+") == 0) {
			direct = 1;
		} else if (strcmp(directText, "-") == 0) {
			direct = -1;
		} else {
			file_error(file, "Invalid value '%s' in 'strand' column", directText);
			file_close(file);
			return -EINVAL;
		}

		if (read_string(file, name, sizeof(name))) {
			file_error(file, "Failed to read 'name' column");
			file_close(file);
			return -EINVAL;
		}
		if (strcmp(ref->map.fragments.data[chrom].name, name) != 0) {
			file_error(file, "Column 'name' does not match");
			file_close(file);
			return -EINVAL;
		}

		if (read_integer(file, &pos)) {
			file_error(file, "Failed to read 'pos' column");
			file_close(file);
			return -EINVAL;
		}
		if (node->pos != pos) {
			file_error(file, "Column 'pos' does not match");
			file_close(file);
			return -EINVAL;
		}

		if (read_integer(file, &size)) {
			file_error(file, "Failed to read 'size' column");
			file_close(file);
			return -EINVAL;
		}
		if (node->size != size) {
			file_error(file, "Column 'size' does not match");
			file_close(file);
			return -EINVAL;
		}

		if (read_integer(file, &uniq)) {
			file_error(file, "Failed to read 'uniq' column");
			file_close(file);
			return -EINVAL;
		}
		ref->index_.data[m].node = ref->nodes.data + index;
		ref->index_.data[m].direct = direct;
		ref->index_.data[m].uniq_count = uniq;
		++m;

		skip_current_line(file);
	}
	assert(m == count);
	ref->index_.size = count;

	file_close(file);
	return 0;
}
