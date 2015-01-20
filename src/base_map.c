#include "base_map.h"

char BASE_MAP[256] = { };
char BASE_REV_MAP[16] = { };

void init_base_map(void)
{
	BASE_MAP['A'] = BASE_MAP['a'] = BASE_A;
	BASE_MAP['C'] = BASE_MAP['c'] = BASE_C;
	BASE_MAP['G'] = BASE_MAP['g'] = BASE_G;
	BASE_MAP['T'] = BASE_MAP['t'] = BASE_T;
	BASE_MAP['U'] = BASE_MAP['u'] = BASE_T;
	BASE_MAP['M'] = BASE_MAP['m'] = BASE_A | BASE_C; /* aMino */
	BASE_MAP['K'] = BASE_MAP['k'] = BASE_G | BASE_T; /* Keto */
	BASE_MAP['R'] = BASE_MAP['r'] = BASE_A | BASE_G; /* puRine */
	BASE_MAP['Y'] = BASE_MAP['y'] = BASE_C | BASE_T; /* pYrimidine */
	BASE_MAP['S'] = BASE_MAP['s'] = BASE_C | BASE_G; /* strong */
	BASE_MAP['W'] = BASE_MAP['w'] = BASE_A | BASE_T; /* weak */
	BASE_MAP['B'] = BASE_MAP['b'] = BASE_C | BASE_G | BASE_T; /* not 'A' */
	BASE_MAP['D'] = BASE_MAP['d'] = BASE_A | BASE_G | BASE_T; /* not 'C' */
	BASE_MAP['H'] = BASE_MAP['h'] = BASE_A | BASE_C | BASE_T; /* not 'G' */
	BASE_MAP['V'] = BASE_MAP['v'] = BASE_A | BASE_C | BASE_G; /* not 'T/U' */
	BASE_MAP['N'] = BASE_MAP['n'] = BASE_N;
	BASE_MAP['X'] = BASE_MAP['x'] = BASE_N;

	BASE_REV_MAP[BASE_A] = 'A';
	BASE_REV_MAP[BASE_C] = 'C';
	BASE_REV_MAP[BASE_G] = 'G';
	BASE_REV_MAP[BASE_T] = 'T';
	BASE_REV_MAP[BASE_A | BASE_C] = 'M';
	BASE_REV_MAP[BASE_G | BASE_T] = 'K';
	BASE_REV_MAP[BASE_A | BASE_G] = 'R';
	BASE_REV_MAP[BASE_C | BASE_T] = 'Y';
	BASE_REV_MAP[BASE_C | BASE_G] = 'S';
	BASE_REV_MAP[BASE_A | BASE_T] = 'W';
	BASE_REV_MAP[BASE_C | BASE_G | BASE_T] = 'B';
	BASE_REV_MAP[BASE_A | BASE_G | BASE_T] = 'D';
	BASE_REV_MAP[BASE_A | BASE_C | BASE_T] = 'H';
	BASE_REV_MAP[BASE_A | BASE_C | BASE_G] = 'V';
	BASE_REV_MAP[BASE_N] = 'N';
}
