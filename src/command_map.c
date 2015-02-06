#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "nick_map.h"
#include "ref_map.h"
#include "bn_file.h"

#define DEF_TOLERANCE 0.1
#define DEF_MIN_MATCH 4

static int verbose = 0;
static double tolerance = DEF_TOLERANCE;
static int min_match = DEF_MIN_MATCH;

static void print_usage(void)
{
	fprintf(stderr, "\n"
			"Usage: bntools map [options] <ref> <query>\n"
			"\n"
			"Options:\n"
			"   <ref>        reference genome, in tsv/cmap format\n"
			"   <query>      query molecules/contigs, in tsv/cmap/bnx format\n"
			"   -e <FLOAT>   tolerance to compare fragment size [%f]\n"
			"   -m <INT>     minimal matched labels in query fragment [%d]\n"
			"   -v           show verbose message\n"
			"   -h           show this help\n"
			"\n", DEF_TOLERANCE, DEF_MIN_MATCH);
}

static void print_header(void)
{
	printf("#name\tchrom\tpos\tstrand\tsize\tlabels\t"
			"qsize\tqlabels\trstart\trend\tqstart\tqend\n");
}

static inline size_t max(size_t a, size_t b) { return (a >= b ? a : b); }

static void output_item(const struct ref_map *ref, const struct fragment *qry,
		size_t rindex, size_t qindex, int direct, size_t matched_nicks)
{
	const struct ref_node *p = &ref->nodes.data[rindex];
	const char *qname = qry->name;
	const char *rname = ref->map.fragments.data[p->chrom].name;
	size_t rstart = (direct > 0 ? rindex : rindex + 1);
	size_t rend = (direct > 0 ? rindex + matched_nicks - 1 : rindex - matched_nicks + 2);
	size_t qstart = qindex;
	size_t qend = qindex + matched_nicks - 1;
	int ref_size = abs(ref->nodes.data[rend].pos - ref->nodes.data[rstart].pos);
	int qry_size = qry->nicks.data[qindex + matched_nicks - 2].pos - qry->nicks.data[qindex - 1].pos;
	int pos = (direct > 0 ? p->pos : (p - matched_nicks + 1)->pos);
	size_t i;

	fprintf(stdout, "%s\t%s\t%d\t%s\t", qname, rname, pos, (direct > 0 ? "+" : "-"));
	fprintf(stdout, "%d\t%zd\t%d\t%zd\t", ref_size, matched_nicks, qry_size, matched_nicks);
	fprintf(stdout, "%zd\t%zd\t%zd\t%zd\t", ref->nodes.data[rstart].label, ref->nodes.data[rend].label, qstart, qend);

	for (i = 0; i + 1 < matched_nicks; ++i) {
		fprintf(stdout, "%s%d:%d", (i == 0 ? "" : "|"),
				ref->nodes.data[rindex + direct * i].size,
				qry->nicks.data[qindex + i].pos - qry->nicks.data[qindex + i - 1].pos);
	}
	fprintf(stdout, "\n");
}

static void map(const struct ref_map *ref, struct fragment *qry_item)
{
	size_t rindex, qindex, i, j;

	if (qry_item->nicks.size < min_match) {
		if (verbose > 1) {
			fprintf(stderr, "Warning: Skip '%s' for less than %d labels!\n", qry_item->name, min_match);
		}
		return;
	}

	for (qindex = 1; qindex < qry_item->nicks.size; ++qindex) {
		const struct nick *p = &qry_item->nicks.data[qindex];
		int fragment_size = p->pos - (p - 1)->pos;
		for (i = 0; i < ref->index_.size; ++i) {
			const struct ref_index *r = &ref->index_.data[i];
			const struct ref_node *n = r->node;
			rindex = n - ref->nodes.data;
			if (n->size < fragment_size * (1 - tolerance)) continue;
			if (n->size > fragment_size * (1 + tolerance)) break;

			for (j = 0; j + 1 < qry_item->nicks.size; ++j) {
				if (j > 0 || verbose > 1) {
					int ref_size = n[j * r->direct].size;
					int qry_size = (p + j)->pos - (p + j - 1)->pos;
					if (j > 0) { /* no need to compare when j == 0, since we started from the matched first interval */
						if (ref_size < qry_size * (1 - tolerance)) break;
						if (ref_size > qry_size * (1 + tolerance)) break;
					}
					if (verbose > 1) {
						fprintf(stderr, "matched interval: rindex = %zd, qindex = %zd, "
								"j = %zd, ref_size = %d, qry_size = %d\n",
								rindex, qindex, j, ref_size, qry_size);
					}
				}
				if (r->direct > 0) {
					if (n->flag & LAST_INTERVAL) {
						++j;
						break;
					}
				} else {
					if (n->flag & FIRST_INTERVAL) {
						++j;
						break;
					}
				}
			}
			if (j + 1 >= min_match) { /* now j equals to matched interval number */
				output_item(ref, qry_item, rindex, qindex, r->direct, j + 1);
			}
		}
	}
}

static int check_options(int argc, char * const argv[])
{
	int c;
	while ((c = getopt(argc, argv, "e:m:avh")) != -1) {
		switch (c) {
		case 'e':
			tolerance = atof(optarg);
			break;
		case 'm':
			min_match = atoi(optarg);
			break;
		case 'v':
			++verbose;
			break;
		case 'h':
			print_usage();
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
	char path[PATH_MAX];
	struct ref_map ref;
	struct nick_map qry;
	size_t i;
	struct stat sb;

	if (check_options(argc, argv)) {
		return 1;
	}
	get_index_filename(argv[optind], path, sizeof(path));

	ref_map_init(&ref);
	if (nick_map_load(&ref.map, argv[optind])) {
		return 1;
	}
	if (stat(path, &sb) == -1 && errno == ENOENT) {
		ref_map_build_index(&ref);
	} else {
		if (ref_map_load(&ref, path)) {
			ref_map_free(&ref);
			return 1;
		}
	}

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
