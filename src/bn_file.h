#ifndef __BN_FILE_H__
#define __BN_FILE_H__

#include "nick_map.h"
#include <zlib.h>

enum file_format {  /* of restriction map */
	FORMAT_UNKNOWN = 0,
	FORMAT_TXT,   /* simple text, with one molecule/contig in each line */
	FORMAT_TSV,   /* tab-separated values, with VCF-like comment header  */
	FORMAT_BNX,   /* .bnx file for molecules, by BioNano inc. */
	FORMAT_CMAP,  /* .cmap file for consensus map, by BioNano inc. */
};

struct bn_file {
	gzFile file;
	int format;       /* enum file_format */
	const char *name; /* filename */
	size_t line;      /* current line */
};

struct bn_file *bn_open(const char *filename);
void bn_close(struct bn_file *fp);

int bn_read_header(struct bn_file *fp, struct nick_map *map);
int bn_read(struct bn_file *fp, struct fragment *f);

int parse_format_text(const char *s);

int nick_map_load(struct nick_map *map, const char *filename);
int nick_map_save(const struct nick_map *map, const char *filename, int format);

struct name_list { array(char *) names; };

int load_name_list(struct name_list *name_list, const char *filename);
void free_name_list(struct name_list *name_list);
int name_list_has(const struct name_list *name_list, const char *name);

gzFile open_gzfile_write(const char *filename);
int save_header(gzFile file, const struct nick_map *map, int format);
int save_fragment(gzFile file, const struct fragment *fragment, int format);

int bn_skip_comment_lines(struct bn_file *fp);
int read_string(struct bn_file *fp, char *buf, size_t bufsize);
int read_integer(struct bn_file *fp, int *value);
void skip_current_line(struct bn_file *fp);

#define bn_file_error(fp, fmt, args...) \
	fprintf(stderr, "Error: " fmt " at line %zd of file '%s'\n", ##args, (fp)->line, (fp)->name)

#endif /* __BN_FILE_H__ */
