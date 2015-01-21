#include <stdio.h>
#include <string.h>
#include "version.h"

extern int nick_main (int argc, char * const argv[]);
extern int view_main (int argc, char * const argv[]);
extern int align_main(int argc, char * const argv[]);
extern int map_main  (int argc, char * const argv[]);

static int version_main(int argc, char * const argv[])
{
	fprintf(stdout, "%s\n", VERSION);
	return 0;
}

struct command {
	const char *name;
	int (*proc)(int argc, char * const argv[]);
	const char *desc;
};

static const struct command CMD[] = {
	{ "version", version_main, "Print program version" },
	{ "nick",    nick_main,    "Generate restriction map from sequence" },
	{ "view",    view_main,    "Convert restriction map between formats" },
	{ "align",   align_main,   "Align between two restriction maps" },
	{ "map",     map_main,     "Map molecules to reference genome" },
};

static void print_usage(void)
{
	size_t i;
	size_t width = 0;
	char fmt[32];

	fprintf(stderr, "\n"
			"Program: bntools (Tools for BioNano data analysis)\n"
			"Version: "VERSION"\n"
			"Author : Linlin Yan (yanll<at>mail.cbi.pku.edu.cn)\n"
			"Copyright: 2014-2015, Centre for Bioinformatics, Peking University, China\n"
			"Website: http://github.com/yanlinlin82/bntools/\n"
			"         http://www.cbi.pku.edu.cn/\n"
			"         http://www.bionanogenomics.com/\n"
			"\n"
			"Usage: bntools <command> [options]\n"
			"\n"
			"Command:\n");

	for (i = 0; i < sizeof(CMD) / sizeof(CMD[0]); ++i) {
		size_t len = strlen(CMD[i].name);
		if (len > width) {
			width = len;
		}
	}
	snprintf(fmt, sizeof(fmt), "   %%-%zds   %%s\n", width);

	for (i = 0; i < sizeof(CMD) / sizeof(CMD[0]); ++i) {
		fprintf(stderr, fmt, CMD[i].name, CMD[i].desc);
	}
	fprintf(stderr, "\n");
}

int main(int argc, char * const argv[])
{
	size_t i;
	if (argc == 1) {
		print_usage();
	} else {
		for (i = 0; i < sizeof(CMD) / sizeof(CMD[0]); ++i) {
			if (strcmp(argv[1], CMD[i].name) == 0) {
				return (*CMD[i].proc)(argc - 1, argv + 1);
			}
		}
		fprintf(stderr, "Error: Unrecognized command '%s'\n", argv[1]);
	}
	return 1;
}
