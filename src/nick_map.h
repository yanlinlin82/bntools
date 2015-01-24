#ifndef __NICK_MAP_H__
#define __NICK_MAP_H__

#include <stdint.h>
#include <string.h>
#include "array.h"

#define MAX_ENZYME_NAME_SIZE 16
#define MAX_REC_SEQ_SIZE 32
#define MAX_CHROM_NAME_SIZE 64

enum nick_strand {
	STRAND_UNKNOWN = 0,
	STRAND_PLUS    = 1,
	STRAND_MINUS   = 2,
	STRAND_BOTH    = (STRAND_PLUS | STRAND_MINUS),
	STRAND_END     = 4,
};

enum file_format {  /* of restriction map */
	FORMAT_UNKNOWN = 0,
	FORMAT_TXT,   /* simple text, with one molecule/contig in each line */
	FORMAT_TSV,   /* tab-separated values, with VCF-like comment header  */
	FORMAT_BNX,   /* .bnx file for molecules, by BioNano inc. */
	FORMAT_CMAP,  /* .cmap file for consensus map, by BioNano inc. */
};

struct nick {
	int pos;
	int strand;
};

struct fragment {  /* molecule, contig or chromosome */
	char *name;
	int size;  /* in bp */
	array(struct nick) nicks;
};

struct nick_map {
	array(struct fragment) fragments;

	char enzyme[MAX_ENZYME_NAME_SIZE];
	char rec_seq[MAX_REC_SEQ_SIZE];

	char rec_bases[MAX_REC_SEQ_SIZE];
	int rec_seq_size;
	int nick_offset;
	int palindrome;
};

/* basic functions */

void nick_map_init(struct nick_map *map);
void nick_map_free(struct nick_map *map);

int nick_map_set_enzyme(struct nick_map *map, const char *enzyme, const char *rec_seq);

struct fragment *nick_map_add_fragment(struct nick_map *map, const char *name);
int nick_map_add_site(struct fragment *f, int pos, int strand);

/* IO functions */

int parse_format_text(const char *s);

int nick_map_load(struct nick_map *map, const char *filename);
int nick_map_save(const struct nick_map *map, const char *filename, int format);

int nick_map_load_fasta(struct nick_map *map, const char *filename, int chrom_only, int verbose);

#endif /* __NICK_MAP_H__ */
