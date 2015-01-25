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
	free(ref->nodes);
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
		for (strand = STRAND_PLUS, matched = 0; strand <= STRAND_MINUS; ++strand) {
			if (matched || seq_match(buf->data + buf->pos - ref->rec_seq_size, ref->rec_bases, ref->rec_seq_size, strand)) {
				int site_pos = base_count - (strand == STRAND_MINUS ? ref->nick_offset : (ref->rec_seq_size - ref->nick_offset));
				if (nick_map_add_site(f, site_pos, strand)) {
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
				ret = nick_map_add_site(f, base_count, STRAND_END);
				if (ret) {
					goto out;
				}
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
		ret = nick_map_add_site(f, base_count, STRAND_END);
		if (ret) {
			goto out;
		}
		f->size = base_count;
		if (verbose > 0) {
			fprintf(stderr, "%d bp\n", base_count);
		}
	}
out:
	gzclose(file);
	return ret;
}

static int sort_by_size(const void *a, const void *b)
{
	const struct node *pa = *(const struct node * const *)a;
	const struct node *pb = *(const struct node * const *)b;
	for (;;) {
		if (pa->size < pb->size) return -1;
		if (pa->size > pb->size) return -1;
		if ((pa->flag & LAST_INTERVAL) && (pb->flag & LAST_INTERVAL)) return 0;
		if (pa->flag & LAST_INTERVAL) return -1;
		if (pb->flag & LAST_INTERVAL) return 1;
		++pa;
		++pb;
	}
	return 0;
}

void generate_ref_nodes(struct ref_map *ref)
{
	size_t i, j, k;

	assert(ref->size == 0);
	assert(ref->nodes == NULL);
	assert(ref->index == NULL);

	ref->size = 0;
	for (i = 0; i < ref->map.fragments.size; ++i) {
		ref->size += ref->map.fragments.data[i].nicks.size - 1;
	}

	ref->nodes = malloc(sizeof(struct node) * ref->size);
	ref->index = malloc(sizeof(struct node *) * ref->size);

	for (i = 0, k = 0; i < ref->map.fragments.size; ++i) {
		const struct fragment *f = &ref->map.fragments.data[i];
		assert(f->nicks.size > 1);
		assert(f->nicks.data[0].pos == 0);
		assert(f->nicks.data[f->nicks.size - 1].pos == f->size);
		for (j = 1; j < f->nicks.size; ++j) {
			ref->nodes[k].chrom = i;
			ref->nodes[k].pos = f->nicks.data[j].pos;
			ref->nodes[k].size = f->nicks.data[j].pos - f->nicks.data[j - 1].pos;
			ref->nodes[k].flag = (j == 1 ? FIRST_INTERVAL : 0) | (j + 1 == f->nicks.size ? LAST_INTERVAL : 0);
			ref->nodes[k].uniq_count = 0;
			ref->index[k] = &ref->nodes[k];
			++k;
		}
	}

	qsort(ref->index, ref->size, sizeof(struct node *), sort_by_size);

	for (i = 0; i + 1 < ref->size; ++i) {
		struct node *a = ref->index[i];
		struct node *b = ref->index[i + 1];
		struct node *end = ref->nodes + ref->size;
		for (j = 0; a + j < end && b + j < end; ++j) {
			if (a[j].size != b[j].size) {
				break;
			}
		}
		if (a->uniq_count < j + 1) {
			a->uniq_count = j + 1;
		}
		if (b->uniq_count < j + 1) {
			b->uniq_count = j + 1;
		}
	}
}

void print_sorted_ref(const struct ref_map *ref)
{
	size_t i, j;

	fprintf(stdout, "#id\tname\tpos\tsize\tflag\tuniq\tseq\n");

	for (i = 0; i < ref->size; ++i) {
		fprintf(stdout, "%zd\t%s\t%d\t%d\t%d\t%d\t", i + 1,
				ref->map.fragments.data[ref->index[i]->chrom].name,
				ref->index[i]->pos, ref->index[i]->size, ref->index[i]->flag,
				ref->index[i]->uniq_count);
		for (j = 0; j < ref->index[i]->uniq_count; ++j) {
			fprintf(stdout, "%s%d", (j == 0 ? "": ","), ref->index[i][j].size);
		}
		fprintf(stdout, "\n");
	}
}
