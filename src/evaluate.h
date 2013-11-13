#ifndef EVALUATE_H
#define EVALUATE_H

#include <QObject>
#include <QString>
#include <vector>

#include "event.h"
#include "polynom.h"
#include "addback.h"

#define LEAFCNT 8
#define CALPREFIX "leaf"
// e.g. leaf0.cal
#define STRETCHFACTOR 1
// factor by which calibrated energy value will be multiplied

#define HISTMAX LEAFCNT
#define HISTCHANNELMAX 65536

#define TIMEHISTMAX LEAFCNT*(LEAFCNT-1)/2

#define MATMAX LEAFCNT*(LEAFCNT-1)/2 + 1
#define MATCHANNELMAX 8192
#define MATSTACKMAX 5000000
#define MATBYTESIZE MATCHANNELMAX*MATCHANNELMAX*sizeof(int)

// prefix for histogram files, e.g. singlehist0
#define HISTPREFIX "singles/singlehist"
//#define HISTPREFIX_BGO "singlehist_bgo"
#define HISTPREFIX_TRIPLE "triples/triplehist"
#define HISTPREFIX_QUAD "quads/quadhist"
#define HISTPREFIX_TIME "timespectra/timehist"

// values for time-cuts
#define TCUT 1
#define TR 4

//walk threshold
#define WALK_THRESH 28
// 0 if no walk-correction for the timespectra and banana-diagrams should be used
// else set any value different from 0
#define WALK_CORRECTION 1

class evaluate : public QObject
{
    Q_OBJECT
public:
    explicit evaluate(QObject *parent = 0);    

    bool addback;
    bool createMatrices;
    bool addCoinc;
    bool debug;
    bool countingRate;

    QString countingRateFile;
    std::vector <double> countingRates;
    std::vector <double> countingRatesError;
    std::vector <double> crTimes;
    long int crHits;
    long int crTimeDiff;
    int crTop;
    int crBottom;


    int cloverNum;
    QString calPrefix;

    int histMax;
    int histChannelMax;
    int leafCnt;
    int stretchFactor;
    int timeHistMax;

    std::vector<int> stats_coinc;   //statistics
    std::vector<t_polynom> cals;    //calibration

    //histograms
    std::vector< std::vector <long int> > hists;
    std::vector< std::vector <long int> > hists_bgo;
    std::vector< std::vector <long int> > hists_triple;
    std::vector< std::vector <long int> > hists_quad;

    //timespectra
    std::vector< std::vector <long int> > timehist;

    //addback histograms
    ab_hist clover[2];

    //Matrices
    std::vector<int> matrix_fds;
    std::vector< std::vector <int> > matrix_stack;
    std::vector<int> matrix_stack_index;
    std::vector<int> matrix_created;
    int32_t maTrix[MATCHANNELMAX][MATCHANNELMAX];
    
signals:
    
public slots:

    int eveInit();
    void get_matrixfilename(char *dest, int mat);
    int get_timespec_index(int detA, int detB);
    int get_matrix_index(int detA, int detB);
    void create_matrixfile(int pos);
    void read_matrix(char *name, int *vals);
    int write_matrix(char *name, int *vals, int rows, int cols);
    void flush_stack(int pos);
    void inc_matrix(int pos, int en_detA, int en_detB);
    void add_coincidence(int detA, int enA, t_polynom calA, int detB, int enB, t_polynom calB, int veto);
    void write_hist(char *name, std::vector<long int> hist, int len);
    void write_hist(char *name, long int* hist, int len);
    void gettestevent(t_event *ev);
    double walkcorrect(double en);
    int sort(int arc, char** agrv);
    
};

#endif // EVALUATE_H
