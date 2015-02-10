#ifndef __BN_FILE_H__
#define __BN_FILE_H__

#include "nick_map.h"
#include "io_base.h"

enum file_format {  /* of restriction map */
	FORMAT_UNKNOWN = 0,
	FORMAT_TXT,   /* simple text, with one molecule/contig in each line */
	FORMAT_TSV,   /* tab-separated values, with VCF-like comment header  */
	FORMAT_BNX,   /* .bnx file for molecules, by BioNano inc. */
	FORMAT_CMAP,  /* .cmap file for consensus map, by BioNano inc. */
};

int bn_read_header(struct file *fp, int *format, struct nick_map *map);
int bn_read(struct file *fp, int format, struct fragment *f);

int parse_format_text(const char *s);

int nick_map_load(struct nick_map *map, const char *filename);
int nick_map_save(const struct nick_map *map, const char *filename, int format);

gzFile open_gzfile_write(const char *filename);
int save_header(gzFile file, const struct nick_map *map, int format);
int save_fragment(gzFile file, const struct fragment *fragment, int format);

int bn_skip_comment_lines(struct file *fp);

#endif /* __BN_FILE_H__ */
