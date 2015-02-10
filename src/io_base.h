#ifndef __IO_BASE_H__
#define __IO_BASE_H__

#include <stdio.h>
#include <zlib.h>

struct file {
	gzFile file;
	const char *name; /* filename */
	size_t line;      /* current line */
};

struct file *file_open(const char *filename);
void file_close(struct file *fp);

static inline int current_char(struct file *fp)
{
	return gzungetc(gzgetc(fp->file), fp->file);
}

void skip_spaces(struct file *fp);
void skip_current_line(struct file *fp);

int read_string(struct file *fp, char *buf, size_t bufsize);
int read_integer(struct file *fp, int *value);

#define file_error(fp, fmt, args...) \
	fprintf(stderr, "Error: " fmt " at line %zd of file '%s'\n", \
			##args, (fp)->line, (fp)->name)

#endif /* __IO_BASE_H__ */
