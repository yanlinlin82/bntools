#include <stdio.h>
#include <string.h>
#include "version.h"
#include "base_map.h"

extern int nick_main (int argc, char * const argv[]);
extern int view_main (int argc, char * const argv[]);
extern int stat_main (int argc, char * const argv[]);
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
	{ "version", version_main, "print program version" },
	{ "nick",    nick_main,    "generate restriction map from sequence" },
	{ "view",    view_main,    "convert restriction map between formats" },
	{ "stat",    stat_main,    "stat motifs in restriction map" },
	{ "align",   align_main,   "align between two restriction maps" },
	{ "map",     map_main,     "map molecules to reference genome" },
};

static void print_usage(void)
{
	size_t i;
	size_t width = 0;
	char fmt[32];

	fprintf(stderr, "\n"
			"Program: bntools (Tools for BioNano data analysis)\n"
			"Version: "VERSION"\n"
			"\n"
			"Usage: bntools <command> [options]\n"
			"\n"
			"Commands:\n");

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
		base_map_init();
		for (i = 0; i < sizeof(CMD) / sizeof(CMD[0]); ++i) {
			if (strcmp(argv[1], CMD[i].name) == 0) {
				return (*CMD[i].proc)(argc - 1, argv + 1);
			}
		}
		fprintf(stderr, "Error: Unrecognized command '%s'\n", argv[1]);
	}
	return 1;
}
