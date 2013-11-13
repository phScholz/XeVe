#ifndef EVENT_H
#define EVENT_H

#include <stdint.h>

#define DETMAX 32
//#define CHANNELMAX 16384
//#define CHANNELMAX 8192
//#define CHANNELMAX 8
#define EVENT_MAXHITS 32

typedef struct {
	int cnt;
	int64_t times[EVENT_MAXHITS];
	int32_t energies[EVENT_MAXHITS];
	int dets[DETMAX];
	int veto[DETMAX];
} t_event;
#endif
