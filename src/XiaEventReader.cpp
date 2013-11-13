#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "dgf_format.h"
#include "event.h"
#include "XiaEventReader.h"
#include "advanced_strings.h"

/*
   written by: Michael Elvers elvers@ikp.uni-koeln.de
   parts of the code adopted from Nigel Warr
*/

XiaEventReader::XiaEventReader() {
	for(int i=0; i<EVENT_TLHISTMAX; i++) {
		evtl_hist[i] = 0;
	}
	sprintf(detectorfile, XER_CONFIGFILE);

	detcnt = 0;
}

XiaEventReader::~XiaEventReader() {
	printf("writing timedifference to histogram\n");
	FILE *f = fopen(EVENT_TLHISTNAME, "w");
	if(!f) {
		printf("cannot open %s\n", EVENT_TLHISTNAME);
			return;
	}

	char info[2000];
	sprintf(info, "subevent coincidence window was set to %d (%lf ns)", EVENT_TIMELIMIT, EVENT_TIMELIMIT_FACTOR*EVENT_TIMELIMIT);
	printf("%s\n", info);


	fprintf(f, "#%s\n", info);
	fprintf(f, "#ticks time/ns counts\n");
	for(int i=0; i<EVENT_TLHISTMAX; i++) {
		fprintf(f, "%d %lf %d\n", i, EVENT_TIMELIMIT_FACTOR*i, evtl_hist[i]);
	}

	fclose(f);
}

// format: detector module subaddress
void XiaEventReader::read_detectorfile(const char *filename) {
	FILE *f = fopen(filename, "r");
	if(!f) {
		fprintf(stderr, "cannot open %s\n", filename);
		exit(1);
	}

	for(int i=0; i<DETMAX; i++) {
		detectors[i].module = 0;
		detectors[i].subaddress = 0;
		detectors[i].nr = 0;
	}

	char buf[2000];
	detcnt = 0;
	while(fgets(buf, 2000, f)) {
		if(as_is_comment(buf)) continue;

		char *tr = as_trim(buf);
		if(strlen(tr)>0) {
			if(detcnt>=DETMAX) {
				fprintf(stderr, "too many detectors, maximum is %d\n", DETMAX);
				exit(1);
			}
			// remove '\n'
			int pos = strlen(tr)-1;
			if(tr[pos]=='\n') {
				tr[pos] = '\0';
			}

			sscanf(tr, "%d %d %d", &detectors[detcnt].nr, &detectors[detcnt].module, &detectors[detcnt].subaddress);
			detcnt++;
		}
	}
	fclose(f);
}

void XiaEventReader::reset() {
	for(int i=0; i<XER_MODMAX; i++) {
		dgf_buffers[i].evcnt = 0;
	}
	currenttime = -1;
	dgf_buffer_index = 0;
}


int XiaEventReader::read_masterheader(FILE *f) {
	uint32_t header[2];
	int pos = ftell(f);
	int firstheader = 1;
	while(1) {
		// read header
		if(fread(header, sizeof(uint32_t), 2, f)!=2) {
			fprintf(stderr, "cannot read master header, file too short\n");
			return 0;
		}

		pos += header[1];
		// extract date information
		if(header[0]!=XER_DATA_HEADER&&firstheader) {
			void *p = malloc(header[1]);
			if(!p) {
                printf("read_masterheader bad alloc");
				fprintf(stderr, "warning: cannot allocate %u bytes. No header information available\n", header[1]);                
			}
			else {
				fread(p, sizeof(char), header[1], f);
				char *str = (char *)p;
				str = &str[42];
				str[strlen(str)-1]='\0';
				strncpy(infostr, str, 2000);
				free(p);
			}
			firstheader = 0;
		}
		// advance to next header
		fseek(f, pos, SEEK_SET);
		
		// check if data header
		if(header[0]==XER_DATA_HEADER) return 1;
	}

	return 0;
}



int64_t XiaEventReader::abs_int64(int64_t val) {
	return (val>0) ? val : -val;
}

uint64_t XiaEventReader::get32bittime(uint16_t hi, uint16_t lo) {
	uint64_t time = hi;
	time = time << 16;
	uint64_t buf = lo;
	time = time | buf;

	return time;
}

uint64_t XiaEventReader::get48bittime(uint16_t hi, uint16_t mi, uint16_t lo) {
	uint64_t time = hi;
	time = time << 16;
	uint64_t buf = mi;
	time = time | buf;
	time = time << 16;
	buf = lo;
	time = time | buf;

	return time;

}

uint64_t XiaEventReader::get_eventtime(int module) {
	t_dgf_buffer *buf = &dgf_buffers[module];
	uint64_t time = buf->events[buf->index].time;
	return time;
}

// return the buffer index with smallest eventtime
int XiaEventReader::get_smallesttime_index() {
	uint64_t st = 0;
	int index = -1;
	for(int i=0; i<XER_MODMAX; i++) {
		t_dgf_buffer *buf = &dgf_buffers[i];
		if(buf->index>=buf->evcnt) continue;
		uint64_t nt = get_eventtime(i);
		if((nt<st)||(index==-1)) {
			st = nt;
			index = i;
		}
	}
	return index;
}

// return the index of buffer which is empty/can be used again
int XiaEventReader::get_emptybuffer() {
	int val = dgf_buffer_index;
	if(dgf_buffer_index>=XER_MODMAX) {
		fprintf(stderr, "unexpected error: dgf_buffer_index >= XER_MODMAX\n");
		exit(1);
	}
	dgf_buffer_index++;

	return val;
}

// calculates the global detector number
int XiaEventReader::getdet(int module, int channel) {
	for(int i=0; i<DETMAX; i++) {
		if((detectors[i].module==module)&&(detectors[i].subaddress==channel)) {
			//return i;
			return detectors[i].nr;
		}
	}
	return -1;
}

void XiaEventReader::make_events() {
	// set all indices to zero
	for(int i=0; i<XER_MODMAX; i++) {
		dgf_buffers[i].index = 0;
	}
	events_index = 0;
	events_cnt = 0;

	int st_index;
	while((st_index = get_smallesttime_index())!=-1) {
		uint64_t st_val = get_eventtime(st_index);


		// reset global event
		t_event *glev = &events[events_cnt];
		t_event_unixtime *glev_ut = &events_unixtime[events_cnt];
		glev_ut->value = 0;
		glev_ut->active = 0;
		glev->cnt = 0;
		for(int mod=0; mod<XER_MODMAX; mod++) {
			t_dgf_buffer *buf = &dgf_buffers[mod];

			// reached end of buffer
			if(buf->index>=buf->evcnt) {
				continue;
			}

			uint64_t val = get_eventtime(mod);
			int64_t diff = val-st_val;
			if(diff>=0&&diff<EVENT_TLHISTMAX) {
				evtl_hist[diff]++;
			}
			t_dgf_event *dgfev = &buf->events[buf->index];

			// timestamps do not match
			if(val-st_val>EVENT_TIMELIMIT) continue;

			// add to global event	
			for(int j=0; j<dgfev->cnt; j++) {
				int ch = dgfev->dets[j];
				uint64_t time = dgfev->chtime[j];
				uint16_t en = dgfev->chen[j];				
				int veto = dgfev->veto[j];


				// add hit
				int det = getdet(buf->module, ch);
				if(det==-1) {
					fprintf(stderr, "warning: not detector specified for module %d, channel %d\n", buf->module, ch);
					continue;
				}
				int cnt = glev->cnt;
				glev->dets[cnt] = det;
				glev->times[cnt] = (int64_t)time; // this should work since time is only 48 bit unsigned
				glev->energies[cnt] = en;
				glev->veto[cnt] = veto;
				glev->cnt = glev->cnt + 1;
			}
			if(dgfev->utimefound) {
				if(glev_ut->active) {
					fprintf(stderr, "warning unixtime already defined\n");
				}	
				glev_ut->active = 1;
				glev_ut->value = dgfev->utime;
			}
			buf->index = buf->index + 1;			
		}

		events_cnt++;
	}
}

void XiaEventReader::reset_current_dgfevent(int bufindex) {
	t_dgf_buffer *buf = &dgf_buffers[bufindex];
	buf->events[buf->evcnt].cnt = 0;
	buf->events[buf->evcnt].utimefound = 0;
}

void XiaEventReader::set_dgfeventtime(int bufindex, uint64_t time) {
	t_dgf_buffer *buf = &dgf_buffers[bufindex];
	if(buf->evcnt>=XER_EVENTMAX) {
		fprintf(stderr, "number of events in buffer %d exceed maximum number (%d)\n", bufindex, XER_EVENTMAX);
		exit(1);
	}
	buf->events[buf->evcnt].time = time;
}

void XiaEventReader::add_hit2dgfevent(int bufindex, int det, uint16_t en, uint64_t time, int veto) {
	t_dgf_buffer *buf = &dgf_buffers[bufindex];
	if(buf->evcnt>=XER_EVENTMAX) {
		fprintf(stderr, "number of events in buffer %d exceed maximum number (%d)\n", bufindex, XER_EVENTMAX);
		exit(1);
	}
	t_dgf_event *ev = &buf->events[buf->evcnt];
	int cnt = ev->cnt;
	ev->dets[cnt] = det;
	ev->chen[cnt] = en;
	ev->chtime[cnt] = time;
	ev->veto[cnt] = veto;
	if(unix_timestamp_found) {
		ev->utime = unix_timestamp;
//printf("time: %u\n", unix_timestamp);
		ev->utimefound = unix_timestamp_found;
	}
	/*ev->utime = unix_timestamp;
	ev->utimefound = unix_timestamp_found;*/
	ev->cnt++;
	unix_timestamp_found = 0;
}

void XiaEventReader::showbin(uint16_t val) {
	uint16_t tmp = 1;
	tmp = tmp << 15;
	for(int i=0; i<16; i++) {
		uint16_t buf = val & tmp;
		if(buf) {
			printf("1");
		}
		else {
			printf("0");
		}
		tmp = tmp >> 1;
	}
	printf("\n");
}

int XiaEventReader::treat_buffer(uint16_t *buffer) {
	//DGF_MAXBUF
	uint16_t nwords = buffer[DGF_H_NWORDS];
	uint16_t module = buffer[DGF_H_MODULE];
	uint16_t fmt = buffer[DGF_H_FORMAT] & DGF_FORMAT_MASK;
	uint16_t timehi = buffer[DGF_H_TIMESTAMP];
	uint16_t timemi = buffer[DGF_H_TIMESTAMP+1];
	uint16_t timelo = buffer[DGF_H_TIMESTAMP+2];


	// calc time
	uint64_t buftime = get48bittime(timehi, timemi, timelo);
	
	if(fmt==DGF_FORMAT_UNIX) {
		uint16_t *evbuf = &buffer[DGF_BUFHEADLEN];
		if(evbuf[0]!=1) {
			fprintf(stderr, "incorrect start of eventbuffer for unixtimestamp, 1 expected but %u found\n", evbuf[0]);
			exit(1);
		}
		uint16_t *utimebuf = &evbuf[DGF_EVENTHEADLEN];
		unix_timestamp = (utimebuf[1] << 16) | utimebuf[0];
		unix_timestamp_found = 1;

		return 1;
	}

	// check if current buffer corresponds to previous buffers
	int64_t diff = buftime - currenttime;
	if((currenttime!=-1)&&(abs_int64(diff)>20)) {
		// make events
		make_events();
		return 0;
	}	
	currenttime = (int64_t)buftime;

	// get index of empty buffer
	int bufindex = get_emptybuffer();

	// set module nr
	dgf_buffers[bufindex].module = module;

	uint16_t *evbuf = &buffer[DGF_BUFHEADLEN];
	while((evbuf-buffer)<nwords) {
		uint16_t pat = evbuf[DGF_E_PATTERN];
		uint16_t evtimemi = evbuf[DGF_E_TIMESTAMP];
		uint16_t evtimelo = evbuf[DGF_E_TIMESTAMP+1];
		evbuf = &evbuf[DGF_EVENTHEADLEN];

		uint64_t evtime = get48bittime(timehi, evtimemi, evtimelo);
		// check if evtime is smaller than buffertime
		if(evtime<buftime) {
			timehi++;
			evtime = get48bittime(timehi, evtimemi, evtimelo);
		}

		reset_current_dgfevent(bufindex);
		set_dgfeventtime(bufindex, evtime);
		for(int i=0; i<DGF_NCHAN; i++) {
			// check if channel i had a hit
			if(!(pat & (1 << i))) {
				//printf("pat %d\n", i);
				continue;
			}

			uint16_t chtimelo=0;
			uint16_t en=0;
			uint16_t ndat=0;
			switch(fmt) {
				// long format
				case 0x0100: case 0x1100: case 0x0101: case 0x1101:
					ndat = evbuf[DGF_C_LNWORDS];
					chtimelo = evbuf[DGF_C_LTRIGGER_TIME];
					en = evbuf[DGF_C_LENERGY];
					evbuf = evbuf + ndat;
				break;

				// PSA format
				case 0x0102: case 0x1102:
					chtimelo = evbuf[DGF_C_PTRIGGER_TIME];
					en = evbuf[DGF_C_PENERGY];
					evbuf = evbuf + DGF_PCHANHEADLEN;
				break;

				// short format
				case 0x0103: case 0x1103:
					chtimelo = evbuf[DGF_C_STRIGGER_TIME];
					en = evbuf[DGF_C_SENERGY];
					evbuf = evbuf + DGF_SCHANHEADLEN;
				break;
				default:
					fprintf(stderr, "unknown format: 0x%04x\n", fmt);
					exit(1);

			}

			// calculate
			uint64_t chtime = get48bittime(timehi, evtimemi, chtimelo);
			if(chtime<evtime) {
				uint64_t buf1 = get32bittime(timehi, evtimemi);
				uint64_t buf2 = chtimelo;
				buf1++;
				chtime = buf1 << 16;
				chtime = chtime | buf2;
			}

			// advance to next event
			evbuf = &evbuf[ndat];

			// add to event array
			int veto = (pat & (1 << (i+12))) ? 1 : 0;
			add_hit2dgfevent(bufindex, i, en, chtime, veto);

		}
		dgf_buffers[bufindex].evcnt++;
	}
	return 1;

}


void XiaEventReader::nextRun(const char *name) {
	if(fbuffer) {
		fclose(fbuffer);
	}

	fbuffer = fopen(name, "r");
	if(!fbuffer) {
		fprintf(stderr, "cannot open %s\n", name);
		exit(1);
	}

	read_detectorfile(detectorfile);
	reset();
	read_masterheader(fbuffer);
	bufferset = 0;
	unix_timestamp = 0;
	unix_timestamp_found = 0;
}

void XiaEventReader::setDetectorFile(const char *name) {
	strcpy(detectorfile, name);
}

bool XiaEventReader::hasMoreEvents() {
	if(events_index<events_cnt) {
		return 1;
	}

	reset();
	bool found = 0;
	if(bufferset) {
		treat_buffer(lastbuffer);
		found = 1;
	}

	bufferset = 1;
	int len = DGF_H_NWORDS + 1;
	int status;
	while((status = fread(lastbuffer, sizeof(uint16_t), len, fbuffer))>0) {
		found = 1;
		// read the DGF_H_NWORDS and find the size of the buffer
		if(status!=len) {
			long pos = ftell(fbuffer);
			fprintf(stderr, "warning: filesize too short for buffer header, position %ld, %d words required but only %d available\n", pos, len, status);
			continue;
			//exit(1);
		}
		uint16_t nwords = lastbuffer[DGF_H_NWORDS];
		int datlen = nwords-len;
		status = fread(&lastbuffer[len], sizeof(uint16_t), datlen, fbuffer);
		if(status!=datlen) {
			long pos = ftell(fbuffer);
			fprintf(stderr, "warning: filesize too short for buffer data, position %ld, %d words required but only %d available\n", pos, datlen, status);
			continue;
			//exit(1);
		}
		// buffer does not have the same time as previous buffers
		if(!treat_buffer(lastbuffer)) {
			return 1;
		}
	}
	bufferset = 0;

	// end of file reached --> make last events
	make_events();

	return found;	
}

int XiaEventReader::readEvent(t_event *dest) {
	t_event_unixtime buf;
	return readEventAndTime(dest, &buf);
}

int XiaEventReader::readEventAndTime(t_event *dest, t_event_unixtime *dest_ut) {
	if(hasMoreEvents()) {
		// copy event
		memcpy(dest, &events[events_index], sizeof(t_event));
		// copy event unixtime
		memcpy(dest_ut, &events_unixtime[events_index], sizeof(t_event_unixtime));
		events_index++;
		return 0;
	}
	else {
		return 1;
	}
}

void XiaEventReader::getInfo(char *dest) {
	//sprintf(dest, "no information implemented yet");
	strcpy(dest, infostr);
}

int XiaEventReader::getDataSkipped() {
	return 0;
}

int XiaEventReader::getDataTotal() {
	int pos = (int)ftell(fbuffer);
	return pos;
}

int XiaEventReader::getDetMax() {
	int max = 0;
	for(int i=0; i<detcnt; i++) {
		if(max<detectors[i].nr) {
			max = detectors[i].nr;
		}
	}
	return max;
}

