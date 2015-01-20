#ifndef __BASE_MAP_H__
#define __BASE_MAP_H__

#define BASE_A 1
#define BASE_C 2
#define BASE_G 4
#define BASE_T 8
#define BASE_N (BASE_A | BASE_C | BASE_G | BASE_T)

extern char BASE_MAP[256];
extern char BASE_REV_MAP[16];

void init_base_map(void);

#endif /* __BASE_MAP_H__ */
