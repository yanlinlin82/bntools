#include <stdio.h>
#include <unistd.h>
#include "nick_map.h"
#include "ref_map.h"

#define DEF_TOLERANCE 0.1
#define DEF_MIN_MATCH 3

static int verbose = 0;
static double tolerance = DEF_TOLERANCE;
static int min_match = DEF_MIN_MATCH;
static int output_align = 0;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools map [options] <ref> <query>\n"
			"\n"
			"Options:\n"
			"   <ref>        reference genome, in tsv/cmap format\n"
			"   <query>      query molecules/contigs, in tsv/cmap/bnx format\n"
			"   -v           show verbose message\n"
			"   -e <FLOAT>   tolerance to compare fragment size [%f]\n"
			"   -m <INT>     minimal matched fragment [%d]\n"
			"   -a           output alignment detail\n"
			"\n", DEF_TOLERANCE, DEF_MIN_MATCH);
}

static void print_header(void)
{
	printf("#name\tchrom\tpos\tstrand\tsize\tlabels\t"
			"qsize\tqlabels\trstart\trend\tqstart\tqend\n");
}

static inline size_t max(size_t a, size_t b) { return (a >= b ? a : b); }

static int get_number_width(int x)
{
	int width = 1;
	while (x >= 10) {
		x /= 10;
		++width;
	}
	return width;
}

static void output_item(const struct ref_map *ref, const struct fragment *qry,
		size_t rindex, size_t qindex, int direct, size_t count)
{
	const char *qname = qry->name;
	const char *rname = ref->map.fragments.data[ref->nodes[rindex].chrom].name;
	size_t rstart = (direct > 0 ? rindex : rindex - count + 1);
	size_t rend = (direct > 0 ? rindex + count : rindex + 1);
	size_t qstart = (direct > 0 ? qindex : qindex + count - 1);
	size_t qend = (direct > 0 ? qindex + count : qindex - 1);

	printf("%s\t", qname);
	printf("%s\t", rname);
	printf("%d\t", ref->nodes[rstart].pos);
	printf("%s\t", (direct > 0 ? "+" : "-"));
	printf("%d\t", ref->nodes[rend].pos - ref->nodes[rstart].pos);
	printf("%zd\t", count);
	printf("%d\t", qry->_nicks.data[qindex + count].pos - qry->_nicks.data[qindex].pos);
	printf("%zd\t", count);
	printf("%zd\t", rstart);
	printf("%zd\t", rend);
	printf("%zd\t", qstart);
	printf("%zd\n", qend);

	if (output_align) {
		int name_width, pos_width;
		size_t i;
		int *q = malloc(sizeof(int) * count);
		int *r = malloc(sizeof(int) * count);
		int *w = malloc(sizeof(int) * count);

		name_width = (int)max(strlen(rname), strlen(qname));
		pos_width = get_number_width(max(max(rstart, rend), max(qstart, qend)));

		for (i = 0; i < count; ++i) {
			r[i] = ref->nodes[rstart + i].pos - ref->nodes[rstart + i - 1].pos;
			q[i] = qry->_nicks.data[qstart + direct * i].pos
					- qry->_nicks.data[qstart + direct * i - 1].pos;
			w[i] = get_number_width(max(r[i], q[i]));
		}

		printf(" %*s: ", name_width, rname);
		for (i = 0; i < count; ++i) {
			printf("%*zd", pos_width, rstart + i);
			if (i + 1 < count) {
				printf(" -(%*d)- ", w[i], r[i]);
			}
		}
		printf("\n");

		printf(" %*s  ", name_width, "");
		for (i = 0; i < count; ++i) {
			printf("%*c", pos_width, '|');
			if (i + 1 < count) {
				printf("  (%*d)  ", w[i], q[i] - r[i]);
			}
		}
		printf("\n");

		printf(" %*s: ", name_width, qname);
		for (i = 0; i < count; ++i) {
			printf("%*zd", pos_width, qstart + direct * i);
			if (i + 1 < count) {
				printf(" -(%*d)- ", w[i], q[i]);
			}
		}
		printf("\n");

		free(w);
		free(r);
		free(q);
	}
}

static void map(const struct ref_map *ref, struct fragment *qry_item)
{
	size_t rindex, qindex, i, j;
	int direct;

	if (qry_item->_nicks.size < min_match) {
		if (verbose > 1) {
			fprintf(stderr, "Warning: Skip '%s' for less than %d labels!\n", qry_item->name, min_match);
		}
		return;
	}

	for (qindex = 1; qindex < qry_item->_nicks.size; ) {
		const struct nick *p = &qry_item->_nicks.data[qindex];
		int fragment_size = p->pos - (p - 1)->pos;
		size_t max_count = 0;
		for (i = 0; i < ref->size; ++i) {
			if (ref->index[i]->size < fragment_size * (1 - tolerance)) continue;
			if (ref->index[i]->size > fragment_size * (1 + tolerance)) break;
			rindex = ref->index[i] - ref->nodes;
			for (direct = 1; direct >= -1; direct -= 2) {
				if (qry_item->_nicks.size > 0) {
					if (direct > 0 && (ref->index[i]->flag | LAST_INTERVAL) != 0) {
						continue;
					} else if (direct < 0 && (ref->index[i]->flag | FIRST_INTERVAL) != 0) {
						continue;
					}
				}
				for (j = 1; j < qry_item->_nicks.size; ++j) {
					int ref_size, qry_size;
					ref_size = ref->nodes[rindex + direct * j].size;
					qry_size = (p + j)->pos - (p + j - 1)->pos;
					if (ref_size < qry_size * (1 - tolerance)) break;
					if (ref_size > qry_size * (1 + tolerance)) break;
					if (direct > 0 && (ref->index[i]->flag | LAST_INTERVAL)) break;
					if (direct < 0 && (ref->index[i]->flag | FIRST_INTERVAL)) break;
				}
				if (j > min_match) {
					output_item(ref, qry_item, rindex, qindex, direct, j);
					if (max_count < j) {
						max_count = j;
					}
				}
			}
		}
		qindex += (max_count > 0 ? max_count : 1);
	}
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "ve:m:a")) != -1) {
		switch (c) {
		case 'v':
			++verbose;
			break;
		case 'e':
			tolerance = atof(optarg);
			break;
		case 'm':
			min_match = atoi(optarg);
			break;
		case 'a':
			output_align = 1;
			break;
		default:
			return 1;
		}
	}
	if (optind + 2 != argc) {
		print_usage();
		return 1;
	}
	return 0;
}

int map_main(int argc, char * const argv[])
{
	struct ref_map ref;
	struct nick_map qry;
	size_t i;

	if (check_options(argc, argv)) {
		return 1;
	}

	ref_map_init(&ref);
	if (nick_map_load(&ref.map, argv[optind])) {
		return 1;
	}
	ref_map_build_index(&ref);

	nick_map_init(&qry);
	if (nick_map_load(&qry, argv[optind + 1])) {
		ref_map_free(&ref);
		return 1;
	}

	print_header();
	for (i = 0; i < qry.fragments.size; ++i) {
		map(&ref, &qry.fragments.data[i]);
	}

	nick_map_free(&qry);
	ref_map_free(&ref);
	return 0;
}
