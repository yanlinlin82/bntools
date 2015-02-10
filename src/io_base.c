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
		fprintf(stderr, "Error: Can not open file '%s'!\n", filename);
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
