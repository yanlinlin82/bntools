#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include <zlib.h>
#include "ref_map.h"
#include "base_map.h"

void ref_map_init(struct ref_map *ref)
{
	memset(ref, 0, sizeof(struct ref_map));
	ref->nick_offset = -1;
}

void ref_map_free(struct ref_map *ref)
{
	free(ref->node);
	free(ref->index);
	nick_map_free(&ref->map);
}

int ref_map_set_enzyme(struct ref_map *ref, const char *enzyme, const char *rec_seq)
{
	char checked_rec_seq[sizeof(ref->map.rec_seq)];
	char rec_bases[sizeof(ref->rec_bases)];
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

	if (ref->nick_offset >= 0) { /* merging maps */
		if (strcmp(ref->map.enzyme, enzyme) != 0 ||
				strcmp(ref->map.rec_seq, checked_rec_seq) != 0) {
			fprintf(stderr, "Error: Merging is supported for only single enzyme\n");
			return 1;
		}
	} else {
		nick_map_set_enzyme(&ref->map, enzyme, checked_rec_seq);

		memcpy(ref->rec_bases, rec_bases, sizeof(ref->rec_bases));
		ref->rec_seq_size = rec_seq_size;
		ref->nick_offset = (nick_offset < 0 ? 0 : nick_offset);
		ref->palindrome = palindrome;
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

static int process_line(struct ref_map *ref, struct fragment *f,
		const char *line, int base_count, struct buffer *buf)
{
	const char *p;
	int strand;
	int matched;
	for (p = line; *p; ++p) {
		if (isspace(*p)) {
			continue;
		}
		if (buf->pos >= sizeof(buf->data)) {
			memcpy(buf->data, buf->data + sizeof(buf->data) - ref->rec_seq_size + 1, ref->rec_seq_size - 1);
			buf->pos = ref->rec_seq_size - 1;
		}
		++base_count;
		buf->data[buf->pos++] = char_to_base(*p);
		if (buf->pos < ref->rec_seq_size) {
			continue;
		}
		for (strand = 0, matched = 0; strand < 2; ++strand) {
			if (matched || seq_match(buf->data + buf->pos - ref->rec_seq_size,
					ref->rec_bases, ref->rec_seq_size, strand)) {
				int site_pos = base_count - (strand == 1 ? ref->nick_offset
						: (ref->rec_seq_size - ref->nick_offset));
				if (nick_map_add_site(f, site_pos,
						(strand == 0 ? NICK_PLUS_STRAND : NICK_MINUS_STRAND))) {
					return -ENOMEM;
				}
				matched = ref->palindrome;
			}
		}
	}
	return base_count;
}

int nick_map_load_fasta(struct ref_map *ref, const char *filename, int chrom_only, int verbose)
{
	gzFile file;
	struct fragment *f = NULL;
	char name[MAX_CHROM_NAME_SIZE] = "";
	struct buffer buf = { };
	int c, ret = 0, base_count = 0;

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
		ret = -EINVAL;
		goto out;
	}
	gzungetc(c, file);

	while (!gzeof(file)) {
		char line[256];
		if (!gzgets(file, line, sizeof(line))) break;
		if (line[0] == '>') {
			if (f) {
				f->size = base_count;
				if (verbose > 0) {
					fprintf(stderr, "%d bp\n", base_count);
				}
			}
			if (chrom_only) {
				int skip = 1;
				const char *p = line + 1;
				if (memcmp(p, "chr", 3) == 0) {
					p += 3;
				}
				if (p[0] >= '1' && p[0] <= '9' && isspace(p[1])) {
					skip = 0;
				} else if ((p[0] == '1' || p[0] == '2') &&
						(p[1] >= '0' && p[1] <= '9') && isspace(p[2])) {
					skip = 0;
				} else if ((p[0] == 'X' || p[1] == 'Y') && isspace(p[1])) {
					skip = 0;
				}
				if (skip) {
					f = NULL;
					continue;
				}
			}
			{
				char *p = line + 1;
				while (*p && !isspace(*p)) ++p;
				*p = '\0';
				snprintf(name, sizeof(name), line + 1);
			}
			if (verbose > 0) {
				fprintf(stderr, "Loading fragment '%s' ... ", name);
			}
			base_count = 0;
			f = nick_map_add_fragment(&ref->map, name);
			if (!f) {
				ret = -ENOMEM;
				goto out;
			}
		} else if (f) {
			int n = process_line(ref, f, line, base_count, &buf);
			if (n < 0) {
				ret = n;
				goto out;
			}
			base_count = n;
		}
	}
	if (f) {
		f->size = base_count;
		if (verbose > 0) {
			fprintf(stderr, "%d bp\n", base_count);
		}
	}
out:
	gzclose(file);
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

void ref_map_build_index(struct ref_map *ref)
{
	size_t i, j, k, m, n;

	assert(ref->node_count == 0);
	assert(ref->node == NULL);
	assert(ref->index_count == 0);
	assert(ref->index == NULL);

	ref->node_count = 0;
	ref->index_count = 0;
	for (i = 0; i < ref->map.fragments.size; ++i) {
		if (ref->map.fragments.data[i].nicks.size > 1) {
			ref->node_count += ref->map.fragments.data[i].nicks.size + 1;
			ref->index_count += (ref->map.fragments.data[i].nicks.size - 1) * 2;
		}
	}

	ref->node = malloc(sizeof(struct ref_node) * ref->node_count);
	ref->index = malloc(sizeof(struct ref_index) * ref->index_count);

	for (i = 0, m = 0, n = 0; i < ref->map.fragments.size; ++i) {
		const struct fragment *f = &ref->map.fragments.data[i];
		if (f->nicks.size <= 1) continue;
		ref->node[m].chrom = i;
		ref->node[m].label = 0;
		ref->node[m].pos = 0;
		ref->node[m].size = f->nicks.data[0].pos;
		ref->node[m].flag = FIRST_INTERVAL;
		++m;
		for (j = 0; j + 1 < f->nicks.size; ++j) {
			ref->node[m].chrom = i;
			ref->node[m].label = j + 1;
			ref->node[m].pos = f->nicks.data[j].pos;
			ref->node[m].size = f->nicks.data[j + 1].pos - f->nicks.data[j].pos;
			ref->node[m].flag = 0;
			for (k = 0; k < 2; ++k) {
				ref->index[n].node = &ref->node[m];
				ref->index[n].direct = (k == 0 ? 1 : -1);
				ref->index[n].uniq_count = 0;
				++n;
			}
			++m;
		}
		ref->node[m].chrom = i;
		ref->node[m].label = f->nicks.size;
		ref->node[m].pos = f->nicks.data[f->nicks.size - 1].pos;
		ref->node[m].size = f->size - f->nicks.data[f->nicks.size - 1].pos;
		ref->node[m].flag = LAST_INTERVAL;
	}

	qsort(ref->index, ref->index_count, sizeof(struct ref_index), sort_by_size);

	for (i = 0; i + 1 < ref->index_count; ++i) {
		struct ref_index *a = &ref->index[i];
		struct ref_index *b = &ref->index[i + 1];
		int x, y, z;
		for (x = 0, y = 0, z = 0; ;
				x += ref->index[i].direct, y += ref->index[i + 1].direct, ++z) {
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
}

void ref_map_dump(const struct ref_map *ref)
{
	size_t i, j;

	fprintf(stdout, "#id\tname\tlabel\tstrand\tpos\tsize\tflag\tuniq\tseq\n");

	for (i = 0; i < ref->index_count; ++i) {
		const struct ref_index *r = &ref->index[i];
		fprintf(stdout, "%zd\t%s\t%zd\t%s\t%d\t%d\t%d\t%d\t", i + 1,
				ref->map.fragments.data[r->node->chrom].name,
				r->node->label, (r->direct > 0 ? "+" : "-"),
				r->node->pos, r->node->size, r->node->flag,
				r->uniq_count);
		for (j = 0; j < r->uniq_count; ++j) {
			const struct ref_node *n = &r->node[j * r->direct];
			fprintf(stdout, "%s%d", (j == 0 ? "": ","), n->size);
			if (meet_last(r, j * r->direct)) break;
		}
		fprintf(stdout, "\n");
	}
}
