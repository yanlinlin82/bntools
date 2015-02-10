#include "io_base.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

struct file *file_open(const char *filename)
{
	struct file *fp;

	assert(filename != NULL);

	fp = malloc(sizeof(struct file));
	if (!fp) {
		return NULL;
	}

	if (strcmp(filename, "-") == 0 || strcmp(filename, "stdin") == 0) {
		fp->file = gzdopen(0, "r"); /* stdin */
	} else {
		fp->file = gzopen(filename, "r");
	}
	if (!fp->file) {
		fprintf(stderr, "Error: Can not open file to read: '%s'\n", filename);
		free(fp);
		return NULL;
	}

	fp->name = filename;
	fp->line = 1;
	return fp;
}

void file_close(struct file *fp)
{
	if (fp) {
		if (fp->file) {
			gzclose(fp->file);
		}
		free(fp);
	}
}

void skip_spaces(struct file *fp)
{
	int c;

	assert(fp != NULL);

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '\n') {
			++fp->line;
		} else if (!isspace(c)) {
			gzungetc(c, fp->file);
			break;
		}
	}
}

void skip_current_line(struct file *fp)
{
	int c;

	assert(fp != NULL);

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '\n') {
			++fp->line;
			break;
		}
	}
}

int read_string(struct file *fp, char *buf, size_t bufsize)
{
	size_t i = 0;
	int c;

	assert(fp != NULL);
	assert(buf != NULL);
	assert(bufsize > 0);

	skip_spaces(fp);
	while ((c = gzgetc(fp->file)) != EOF) {
		if (isspace(c)) {
			gzungetc(c, fp->file);
			break;
		} else if (i + 1 < bufsize) {
			buf[i++] = (char)c;
		}
	}
	buf[i] = '\0';
	return (i > 0 ? 0 : -1);
}

int read_integer(struct file *fp, int *value)
{
	int c;

	assert(fp != NULL);
	assert(value != NULL);

	skip_spaces(fp);
	c = gzgetc(fp->file);
	if (!isdigit(c)) {
		gzungetc(c, fp->file);
		return -1;
	}
	*value = c - '0';
	while ((c = gzgetc(fp->file)) != EOF) {
		if (!isdigit(c)) {
			gzungetc(c, fp->file);
			break;
		}
		*value = *value * 10 + (c - '0');
	}
	return 0;
}
