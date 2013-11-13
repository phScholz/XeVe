#ifndef ADDBACK_H
#define ADDBACK_H

#include <vector>
#include "polynom.h"
#include "event.h"

#define AB_WIDTH 50
#define AB_HISTMAX 65536
#define AB_MINRANGE 100

#ifndef HISTCHANNELMAX
#define HISTCHANNELMAX 65536
#endif

#define HISTPREFIX_ADD1 "ab/abhistv1"
#define HISTPREFIX_ADD2 "ab/abhistv2"
#define HISTPREFIX_ADD3 "ab/abhistv3"
#define HISTPREFIX_TOTAL1 "ab/totalv1"
#define HISTPREFIX_TOTAL2 "ab/totalv2"
#define HISTPREFIX_TOTAL3 "ab/totalv3"
#define HISTPREFIX_DIRECT1 "ab/directv1"
#define HISTPREFIX_DIRECT2 "ab/directv2"
#define HISTPREFIX_DIRECT3 "ab/directv3"
#define HISTPREFIX_SINGLES1 "ab/singlesv1"
#define HISTPREFIX_SINGLES2 "ab/singlesv2"
#define HISTPREFIX_SINGLES3 "ab/singlesv3"

#define LIMIT 300

typedef struct {
    long int total[3][HISTCHANNELMAX]; //Totalmode-Spektrum
    long int addback[3][HISTCHANNELMAX]; //Addbackmode-Spektrum
    long int direct[3][HISTCHANNELMAX]; //Directmode-Spektrum
    long int singles[3][HISTCHANNELMAX]; //Singlesmode-Spektrum
    long int stat_addback[3][10]; //Statistikarray
} ab_hist;

void addbackv1(ab_hist *clover, t_event ev, std::vector<t_polynom> cals, int clonr, int begin);
void addbackv2(ab_hist *clover, t_event ev, std::vector<t_polynom> cals, int clonr, int begin);
void addbackv3(ab_hist *clover, t_event ev, std::vector<t_polynom> cals, int clonr, int begin);
void ab_hist_init(ab_hist *clover, int anzahl);
long int ab_calibrate(t_polynom calA, t_polynom calB, long int en);
long int ab_getmax(long int *hist, long int *sum, int count);
void ab_add2(long int *stats, long int *hist, long int *hist2, long int *sum, int x, int y, int z);
void ab_add3(long int *stats, long int *hist, long int *hist2, long int *sum, int a, int b, int c, int d, int e, int f,int g);
int ab_multi(int cnt, int *dets, int begin);

#endif

