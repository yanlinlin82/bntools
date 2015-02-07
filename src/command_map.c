#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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
			"qsize\tqlabels\trstart\trend\tqstart\tqend\tmissing\textra\talignment\n");
}

static inline size_t max(size_t a, size_t b) { return (a >= b ? a : b); }

static void output_item(const struct ref_map *ref, const struct fragment *qry,
		size_t rindex, size_t qindex, int direct, size_t rlabel, size_t qlabel,
		const int *matches, size_t match_count, size_t missing, size_t extra)
{
	const struct ref_node *p = &ref->nodes.data[rindex];
	const char *qname = qry->name;
	const char *rname = ref->map.fragments.data[p->chrom].name;
	size_t rstart = (direct > 0 ? rindex : rindex + 1);
	size_t rend = (direct > 0 ? rindex + rlabel - 1 : rindex - rlabel + 2);
	size_t qstart = qindex;
	size_t qend = qindex + qlabel - 1;
	int ref_size = abs(ref->nodes.data[rend].pos - ref->nodes.data[rstart].pos);
	int qry_size = qry->nicks.data[qindex + qlabel - 2].pos - qry->nicks.data[qindex - 1].pos;
	int pos = (direct > 0 ? p->pos : (p - rlabel + 1)->pos);
	size_t i, j, k;

	fprintf(stdout, "%s\t%s\t%d\t%s\t", qname, rname, pos, (direct > 0 ? "+" : "-"));
	fprintf(stdout, "%d\t%zd\t%d\t%zd\t", ref_size, rlabel, qry_size, qlabel);
	fprintf(stdout, "%zd\t%zd\t%zd\t%zd\t", ref->nodes.data[rstart].label, ref->nodes.data[rend].label, qstart, qend);
	fprintf(stdout, "%zd\t%zd\t", missing, extra);

	for (i = 0, j = 0, k = 0; i < match_count; ++i) {
		if (i > 0) {
			fprintf(stdout, "|");
		}

		fprintf(stdout, "%d", ref->nodes.data[rindex + direct * j++].size);
		if (matches[i] == 2) {
			fprintf(stdout, "+%d", ref->nodes.data[rindex + direct * j++].size);
		}

		fprintf(stdout, ":%d", qry->nicks.data[qindex + k].pos - qry->nicks.data[qindex + k - 1].pos);
		++k;
		if (matches[i] == 3) {
			fprintf(stdout, "+%d", qry->nicks.data[qindex + k].pos - qry->nicks.data[qindex + k - 1].pos);
			++k;
		}
	}
	fprintf(stdout, "\n");
}

static void map(const struct ref_map *ref, struct fragment *qry_item)
{
	size_t rindex, qindex, i, j, k, missing, extra;
	array(int) matches = { };

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

			matches.size = 0;
			for (j = 0, k = 0, missing = 0, extra = 0; j + 1 < qry_item->nicks.size; ++j, ++k) {
				int match = 0;
				int ref_size, qry_size;

				if (r->direct > 0) {
					if (n[j * r->direct].flag & LAST_INTERVAL) {
						assert(j > 0);
						break;
					}
				} else {
					if (n[j * r->direct].flag & FIRST_INTERVAL) {
						assert(j > 0);
						break;
					}
				}
				if (j == 0) {
					match = 1; /* the first interval is always matched */
					if (verbose > 1) {
						ref_size = n[j * r->direct].size;
						qry_size = (p + k)->pos - (p + k - 1)->pos;
					}
				} else {
					/* try matching */
					ref_size = n[j * r->direct].size;
					qry_size = (p + k)->pos - (p + k - 1)->pos;
					if (qry_size >= ref_size * (1 - tolerance) && qry_size <= ref_size * (1 + tolerance)) {
						match = 1;
					}

					if (!match) {
						/* try matching with missing nick */
						ref_size = n[(j + 1) * r->direct].size + n[j * r->direct].size;
						qry_size = (p + k)->pos - (p + k - 1)->pos;
						if (qry_size >= ref_size * (1 - tolerance) && qry_size <= ref_size * (1 + tolerance)) {
							match = 2;
							++missing;
						}
					}

					if (!match) {
						/* try matching with extra nick */
						ref_size = n[j * r->direct].size;
						qry_size = (p + k + 1)->pos - (p + k - 1)->pos;
						if (qry_size >= ref_size * (1 - tolerance) && qry_size <= ref_size * (1 + tolerance)) {
							match = 3;
							++extra;
						}
					}
				}
				if (!match) break;

				if (array_reserve(matches, matches.size + 1)) {
					fprintf(stderr, "Error: Failed to allocate memory!\n");
					array_free(matches);
					return;
				}
				matches.data[matches.size++] = match;

				if (verbose > 1) {
					fprintf(stderr, "matched interval: rindex = %zd, qindex = %zd, match = %d, "
							"j = %zd, k = %zd, ref_size = %d, qry_size = %d\n",
							rindex, qindex, match, j, k, ref_size, qry_size);
				}

				if (match == 2) {
					++j;
				} else if (match == 3) {
					++k;
				}
			}
			if (matches.size + 1 >= min_match) {
				output_item(ref, qry_item, rindex, qindex, r->direct,
						j + 1, k + 1, matches.data, matches.size, missing, extra);
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
