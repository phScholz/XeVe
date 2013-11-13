#ifndef XIAEVENTREADER_H
#define XIAEVENTREADER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "event.h"
#include "dgf_format.h"

#define XER_MASTER_HEADER 0x696E694D
#define XER_DATA_HEADER 0x61746144

#define XER_MODMAX 10
#define XER_EVENTMAX 10000

#define XER_CONFIGFILE "xia.config"

#define EVENT_TIMELIMIT 100
#define EVENT_TIMELIMIT_FACTOR 12.5 // Factor to multiply EVENT_TIMELIMIT to get coincidence window in ns
#define EVENT_TLHISTMAX 10000 //time limit hist max
#define EVENT_TLHISTNAME "timediffhist" // name of the histogram file

typedef struct {
	uint16_t chen[4];
	uint64_t chtime[4];
	int dets[4];
	int veto[4];
	int cnt;
	uint64_t time; // dgftime
	uint32_t utime; // unixtime
	int utimefound; // 1 if a valid time is in utime
} t_dgf_event;

typedef struct {
	t_dgf_event events[XER_EVENTMAX];
	int module; // module nr
	int evcnt;
	int index; // used in make_events
} t_dgf_buffer;

typedef struct {
	int module;
	int subaddress;
	int nr;
} t_detector;

typedef struct {
	uint32_t value;
	int active; // true if a timestamp was received for this event
} t_event_unixtime;




class XiaEventReader {
	public:
		XiaEventReader();
		~XiaEventReader();

		void nextRun(const char *name);
		void setDetectorFile(const char *name);

		bool hasMoreEvents();
		int readEvent(t_event *event);
		int readEventAndTime(t_event *event, t_event_unixtime *event_ut);
		void getInfo(char *dest);
		int getDataSkipped();
		int getDataTotal();
		int getDetMax();
		
	
	private:

		void read_detectorfile(const char *filename);
		void reset();
		int read_masterheader(FILE *f);
		uint64_t get32bittime(uint16_t hi, uint16_t lo);
		uint64_t get48bittime(uint16_t hi, uint16_t mi, uint16_t lo);
		uint64_t get_eventtime(int module);
		int get_smallesttime_index();
		int get_emptybuffer();
		int getdet(int module, int channel);
		void make_events();
		void reset_current_dgfevent(int bufindex);
		void set_dgfeventtime(int bufindex, uint64_t time);
		void add_hit2dgfevent(int bufindex, int det, uint16_t en, uint64_t time, int veto);
		int treat_buffer(uint16_t *buffer);
		
		int64_t abs_int64(int64_t val);
		void showbin(uint16_t val);

		t_detector detectors[DETMAX];

		t_dgf_buffer dgf_buffers[XER_MODMAX];
		int dgf_buffer_index;

		// contains matched events between buffer
		t_event events[XER_EVENTMAX];
		int events_cnt;
		int events_index;
		int64_t currenttime;

		// contains unix time stamp for events above
		t_event_unixtime events_unixtime[XER_EVENTMAX];

		uint16_t lastbuffer[DGF_MAXBUF];
		int bufferset;
		FILE *fbuffer;


		// event time difference
		int evtl_hist[EVENT_TLHISTMAX];

		// unix timestamp
		int unix_timestamp_found;
		uint32_t unix_timestamp;

		char detectorfile[2000];
		char infostr[2000];

		int detcnt;

		


};
#endif
