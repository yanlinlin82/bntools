#ifndef __BASE_MAP_H__
#define __BASE_MAP_H__

#define BASE_A 1
#define BASE_C 2
#define BASE_G 4
#define BASE_T 8
#define BASE_N (BASE_A | BASE_C | BASE_G | BASE_T)

extern char CHAR_TO_BASE[256];
extern char BASE_TO_CHAR[16];
extern char BASE_TO_COMP[16];

extern void init_base_map(void);

static inline char char_to_base(char c) { return CHAR_TO_BASE[(int)c]; }
static inline char base_to_char(char c) { return BASE_TO_CHAR[(int)c]; }
static inline char base_to_comp(char c) { return BASE_TO_COMP[(int)c]; }

#endif /* __BASE_MAP_H__ */
