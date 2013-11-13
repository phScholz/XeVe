/*
 * written by  Michael Elvers elvers@ikp.uni-koeln.de
 * modified/remixed by Philipp Scholz pscholz@ikp.uni-koeln.de
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string>
#include <error.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <QDebug>
#include <QFile>


#include "XiaEventReader.h"
#include "polynom.h"
#include "mfile.h"
#include "addback.h"


#include "evaluate.h"

evaluate::evaluate(QObject *parent) :
    QObject(parent)
{
    addback=false;
    createMatrices=false;
    addCoinc=false;
    cloverNum=2;
    leafCnt=8;
    histMax=leafCnt;
    calPrefix="leaf";
    stretchFactor=1;
    histChannelMax=65536;
    timeHistMax=leafCnt*(leafCnt-1)/2;    
}



// returns filename for temporary matrices
void evaluate::get_matrixfilename(char *dest, int mat) {
    sprintf(dest, "mat%02d.tmp", mat);
}


int evaluate::eveInit(){

    stats_coinc.assign(DETMAX,0);

    for(int i=0; i<leafCnt; i++){
        t_polynom dummy;
        cals.push_back(dummy);
    }

    for(long int i=0; i<histMax; i++){
        std::vector<long int> dummy;
        hists.push_back(dummy);
        hists_triple.push_back(dummy);
        hists_quad.push_back(dummy);
        hists_bgo.push_back(dummy);
        for(long int j=0; j<histChannelMax; j++){
            hists.at(i).push_back(0);
            hists_triple.at(i).push_back(0);
            hists_quad.at(i).push_back(0);
            hists_bgo.at(i).push_back(0);
        }
    }

    for(long int i=0; i<timeHistMax; i++){
        std::vector<long int> dummy;
        timehist.push_back(dummy);
        for(long int j=0; j<histChannelMax; j++){
            timehist.at(i).push_back(0);
        }
    }

    if(createMatrices){
        for(int i=0; i<MATMAX; i++){
            matrix_fds.push_back(0);
            matrix_stack_index.push_back(0);
            matrix_created.push_back(0);
            std::vector<int> dummy;
            matrix_stack.push_back(dummy);
            for(int j=0; j<MATSTACKMAX; j++){
                matrix_stack.at(i).push_back(0);
            }
        }
    }

    if(addback)
        ab_hist_init(clover, cloverNum);

    // read calibrations
    for(int i=0; i<LEAFCNT; i++) {
        char name[2000];
        sprintf(name, "%s%d.cal", CALPREFIX, i);
        cals[i] = read_polynomfile(name);
    }

    for(int i=0; i<LEAFCNT; i++) {
        //pol_print(cals[i]);
    }

    if(access("singles", F_OK)) {
        if(mkdir("singles", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
            fprintf(stderr, "cannot created directory singles");
            exit(1);
        }
    }

    if(access("triples", F_OK)) {
        if(mkdir("triples", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
            fprintf(stderr, "cannot created directory triples");
            exit(1);
        }
    }

    if(access("triples", F_OK)) {
        if(mkdir("triples", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
            fprintf(stderr, "cannot created directory triples");
            exit(1);
        }
    }

    if(access("quads", F_OK)) {
        if(mkdir("quads", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
            fprintf(stderr, "cannot created directory quads");
            exit(1);
        }
    }


    if(createMatrices){
        // create direcetory for matrices if needed
        if(access("matrices", F_OK)) {
            if(mkdir("matrices", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
                fprintf(stderr, "cannot created directory matrices");
                exit(1);
            }
        }

        // create directory for timespectra if needed
        if(access("timespectra", F_OK)) {
            if(mkdir("timespectra", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
                fprintf(stderr, "cannot created directory matrices");
                exit(1);
            }
        }

        printf("allocating matrix stack memory\n");
        // unfortunately this method is very slow, however it's easy to understand
        /*for(int i=0; i<MATMAX; i++) {
            for(int j=0; j<MATCHANNELMAX; j++) {
                for(int k=0; k<MATCHANNELMAX; k++) {
                    matrices[i][j][k] = 0;
                }
            }
        }*/

        for(int mat=0; mat<MATMAX; mat++) {
                    if(mat%5==0) {
                            printf(" %02d ", mat);
                    }
                    printf("."); fflush(stdout);

                    // close if any files are open
                    if(matrix_fds[mat]!=-1) {
                            close(matrix_fds[mat]);
                    }
                    matrix_fds[mat] = -1;
                    matrix_created[mat] = 0;

                    matrix_stack_index[mat] = 0;
        }


        printf("done\n");
    }

    return 0;
}

int evaluate::get_timespec_index(int detA, int detB) {
    if(detA>detB) {
        fprintf(stderr, "error in line %d: detA (%d) > detB (%d)\n", __LINE__, detA, detB);
        exit(1);
    }
    return (detA*LEAFCNT - detA*(detA+1)/2 + detB - detA - 1);
}

int evaluate::get_matrix_index(int detA, int detB) {
    return get_timespec_index(detA, detB) + 1;
}


// creates a temporary matrix
void evaluate::create_matrixfile(int pos) {
        char cname[2000];
        get_matrixfilename(cname, pos);

        // open matrix
        matrix_fds[pos] = open(cname, O_RDWR | O_CREAT | O_TRUNC, 0777);
        if(matrix_fds[pos]<0) {
                fprintf(stderr, "cannot open %s\n", cname);
                exit(1);
        }
}

// reads from a matrix file if it was created
// return an empty matrix otherwise
void evaluate::read_matrix(char *name, int *vals) {
        // open
        MFILE *mat = mopen(name, (char *)"r");
        if(!mat) {
                fprintf(stderr, "cannot open %s\n", name);
                exit(1);
        }

    // read info
    minfo info;
    mgetinfo(mat, &info);
    if(info.lines!=MATCHANNELMAX) {
        fprintf(stderr, "incorrect number of lines in matrix %s, %d expected but %d found\n", name, MATCHANNELMAX, info.lines);
        exit(1);
    }
    if(info.columns!=MATCHANNELMAX) {
        fprintf(stderr, "incorrect number of columns in matrix %s, %d expected but %d found\n", name, MATCHANNELMAX, info.columns);
        exit(1);
    }

    // read
    for(int lin=0; lin<MATCHANNELMAX; lin++) {
        mgetint(mat, &vals[lin*MATCHANNELMAX], 0, lin, 0, MATCHANNELMAX);
    }

    // close
    mclose(mat);
}

// write final matrix
int evaluate::write_matrix(char *name, int *vals, int rows, int cols) {
        // open
        MFILE *mat = mopen(name, (char *)"w");
        if(!mat) {
                fprintf(stderr, "cannot open %s\n", name);
                return 0;
        }

        // set info and format
        minfo info;
        mgetinfo(mat, &info);
        info.levels = 1;
        info.lines = rows;
        info.columns = cols;
        if(msetinfo(mat, &info)!=0) {
                fprintf(stderr, "cannot set info for matrix file %s\n", name);
                return 0;
        }

        msetfmt(mat, (char *)"8k.8k:lc");

        // write
        for(int lin=0; lin<rows; lin++) {
                mputint(mat, &vals[lin*cols], 0, lin, 0, cols);
        }

        // close
        mclose(mat);

        return 1;
}


// writes current stack to temporary matrix file
void evaluate::flush_stack(int pos) {
        int index = matrix_stack_index[pos];
        //int fd = matrix_fds[pos];

        //printf("flushing pos %d, fd %d\n", pos, fd);

    // allocate disk space for matrix
        void *buf = malloc(MATBYTESIZE);
        if(!buf) {
                fprintf(stderr, "cannot allocate memory\n");
                exit(1);
        }


        int32_t *mat;
        mat = (int32_t *)buf;
    char name[2000];
    get_matrixfilename(name, pos);
    if(!matrix_created[pos]) {
        for(int i=0; i<MATCHANNELMAX*MATCHANNELMAX; i++) {
            mat[i] = 0;
        }
    }
    else {
        read_matrix(name, mat);
    }

    // fill matrix
        for(int i=0; i<index; i++) {
                mat[matrix_stack[pos][i]]++;
        }

    // write back to disk
    write_matrix(name, mat, MATCHANNELMAX, MATCHANNELMAX);
    matrix_created[pos]=1;

    // free memory
        free(buf);

    // reset stack
        matrix_stack_index[pos] = 0;

}


// inc the element after checking if it is a valid event
void evaluate::inc_matrix(int pos, int en_detA, int en_detB) {
    if(en_detA<0||en_detA>=MATCHANNELMAX||en_detB<0||en_detB>=MATCHANNELMAX) {
        //fprintf(stderr, "value pair out of range: index %d, enA %d, enB %d\n", index, enA, enB);
        return;
    }

    //matrix_stack[pos][matrix_stack_index[pos]]=en_detB*MATCHANNELMAX+en_detA;
    matrix_stack[pos][matrix_stack_index[pos]]=en_detA*MATCHANNELMAX+en_detB;
    matrix_stack_index[pos]++;
    // flush if stack is full
        if(matrix_stack_index[pos]==MATSTACKMAX) {
                flush_stack(pos);
        }

}

// add a coincidence to matrix
void evaluate::add_coincidence(int detA, int enA, t_polynom calA, int detB, int enB, t_polynom calB, int veto) {

    // discriminate from noise
    if(enA<100||enB<100) {
        return;
    }

    // calibration
    double x = 1.*enA;
    double p = 1.*random()/RAND_MAX;
    x = x + p;
    x = calA.coeffs[0] + calA.coeffs[1]*x + calA.coeffs[2]*x*x + calA.coeffs[3]*x*x*x;
    // stretch
    x = x * STRETCHFACTOR;
    enA = (int)x;

    x = 1.*enB;
    p = 1.*random()/RAND_MAX;
    x = x + p;
    x = calB.coeffs[0] + calB.coeffs[1]*x + calB.coeffs[2]*x*x + calB.coeffs[3]*x*x*x;
    // stretch
    x = x * STRETCHFACTOR;
    enB = (int)x;

    // add to symmetric matrix
    inc_matrix(0, enA, enB);
    inc_matrix(0, enB, enA);
    // determine index
    int index = get_matrix_index(detA, detB);
    // 1--3
    /*if(detA==0) {
        index = detB;
    }
    else if(detA==1) {
        index = detB + 2;
    }
    else if(detA==2) {
        index = 6;
    }
    else {
        fprintf(stderr, "unexpected value pair for detectors: detA %d, detB %d\n", detA, detB);
        return;
    }*/

    inc_matrix(index, enA, enB);
/*
    if (!veto) {
        inc_matrix(index+7, enA, enB);
    }
*/
}


// write history file
void evaluate::write_hist(char *name, std::vector<long int> hist, int len) {
    // open histogram file
    FILE *f = fopen(name, "w");
    if(!f) {
        fprintf(stderr, "cannot open histogram file %s. Not writing to this file\n", name);
        return;
    }

    // write data
    for(int i=0; i<len; i++) {
        fprintf(f, "%ld\n", hist[i]);
    }

    // close
    fclose(f);
}
void evaluate::write_hist(char *name, long int *hist, int len) {
    // open histogram file
    FILE *f = fopen(name, "w");
    if(!f) {
        fprintf(stderr, "cannot open histogram file %s. Not writing to this file\n", name);
        return;
    }

    // write data
    for(int i=0; i<len; i++) {
        fprintf(f, "%ld\n", hist[i]);
    }

    // close
    fclose(f);
}

// generate a dummy event for testing
void evaluate::gettestevent(t_event *ev) {
    static int cnt = 0;
    /*if(cnt===0)*/ {
        ev->cnt=2;
        ev->dets[0] = 0;
        ev->dets[1] = 1;
        ev->energies[0] = 10;
        ev->energies[1] = 20;
    }
    cnt++;
}

// walk-correction function
double evaluate::walkcorrect(double en) {

    // fit function: E(dt) = exp(a(dt)^2+b(dt)+c)
    double a = 0.0398;
    double b = -1.660;
    double c = 20.529;

    if (a<=0) {
        fprintf(stderr, "Negative or Zero argument of 1/sqrt(a). No correction will be applyed\n");
        return 0;
    }
    double x = 1/sqrt(a);

    if (en<=0) {
        fprintf(stderr, "Negative or Zero argument of log(en). No correction will be applyed\t%f\n", en);
        return 0;
    }

    double y = (pow(b,2.0)/(4.0*a))-c;

    if ((log(en)+y)<0) {
        fprintf(stderr, "Negative or Zero argument of sqrt(log(en)+y). No correction will be applyed\n");
        return 0;
    }

    double corr = x*sqrt(log(en)+y);
    return corr;
}


int evaluate::sort(int argc, char *argv[]) {

    eveInit();


    // set seed for random generator to get reproducible results
    srandom(0);


    // check if enough arguments are provided
    if(argc<2) {
        fprintf(stderr, "usage: evaluate <file1> [<file2>]\n");
        fprintf(stderr, "\tdetector configuration: -c <configfile>\n");
        fprintf(stderr, "\tprint index summary: -s\n");
        exit(1);
    }

    char configname[2000];
    sprintf(configname, "config");

    int write_summary = 0;
    char c;
    while((c = getopt(argc, argv, "sc:")) != -1) {
        switch(c) {
            case 'c':
                strcpy(configname, optarg);
            break;
            case 's':
                write_summary = 1;
            break;
            default:
                if(optopt=='c') {
                    fprintf(stderr, "c requires a filename as an argument\n");
                    exit(1);
                }
            break;
        }
    }
    printf("using configfile %s\n", configname);


    XiaEventReader *er = new XiaEventReader();
    if(!er){
        qDebug() << "Bad alloc at XiaEventReader";
    }
    t_event ev;
    t_event_unixtime evut;
    int count = 0;
    double total_time=0;

    for(int i=optind; i<argc; i++) {
        printf("reading %s\n", argv[i]);
        er->setDetectorFile(configname);
        er->nextRun(argv[i]);
        int64_t firsttime = -1;
        int64_t lasttime = -1;


        while(er->hasMoreEvents()) {
            er->readEventAndTime(&ev, &evut);
            lasttime = ev.times[0];
            if(firsttime==-1) {
                firsttime = lasttime;
            }


            if(addCoinc){
                int count=0;
                bool c1=false, c2=false;
                long int erg1=0, erg2=0;

                for(int i=0; i<ev.cnt; i++){
                    if(ev.dets[i]>=0 && ev.dets[i]<4){
                        erg1+=ab_calibrate(cals[ev.dets[i]],cals[0],ev.energies[i]);
                        c1=true;
                    }
                    if(ev.dets[i]>=4 && ev.dets[i]<8){
                        erg2+=ab_calibrate(cals[ev.dets[i]],cals[4],ev.energies[i]);
                        c2=true;
                    }
                }

                if(c2 && c1){
                    ev.dets[0]=0;
                    ev.dets[1]=4;
                    ev.energies[0]=erg1;
                    ev.energies[1]=erg2;
                    ev.cnt=2;
                }
                if(c2 && !c1){
                    ev.dets[0]=4;
                    ev.energies[0]=erg2;
                    ev.cnt=1;
                }
                if(!c2 && c1){
                    ev.dets[0]=0;
                    ev.energies[0]=erg1;
                    ev.cnt=1;
                }

                if(!c1 && !c2) ev.cnt=0;
            }

            if(countingRate){
                long int en=0;
                for(int i=0; i<ev.cnt; i++){
                    en = ev.energies[i];
                    t_polynom cal=cals[ev.dets[i]];

                    if(en > 0){
                        double x=1.*en;
                        double p=1.*random()/RAND_MAX;
                        x=cal.coeffs[0] + cal.coeffs[1]*en;
                        x= x + p;
                        en=(int)x;
                    }
                    else{
                        en=0;
                    }
                }

                if(en>=crBottom && en<=crTop){
                    crHits++;
                }
            }


            // statistics
            stats_coinc[ev.cnt]++;

            count++;
            if(count%100000==0) {
                printf("."); fflush(stdout);
            }

            /*
            // event has timestamp
            if(evut.active) {
                status = abs(status - 1);
                time_t buf = evut.value;
                printf("time %u %s\n",evut.value, asctime(localtime(&buf)));
            }

            if(status==0) continue;
            */

            // histories
            for(int i=0; i<ev.cnt; i++) {
                int det = ev.dets[i];
                int en = ev.energies[i];
                //int veto = ev.veto[i];

                if(det<0||det>=HISTMAX) {
                    fprintf(stderr, "invalid detector nr: %d\n", det);
                    continue;
                }
                if(en<0||en>=HISTCHANNELMAX) {
                    fprintf(stderr, "invalid energie value for det %d: %d\n", det, en);
                    continue;
                }

                hists[det][en]++;
/*
                if (!veto) {
                    hists_bgo[det][en]++;
                }
*/
                // triple coinc hist (EXACT 3!!!)
                if(ev.cnt==3) {
                    hists_triple[det][en]++;
                }

                // quad coinc
                if(ev.cnt==4) {
                    hists_quad[det][en]++;
                }
            }

            if(addback){
                //old ADDBACK
                addbackv1(clover, ev, cals, 0, 0);
                addbackv1(clover, ev, cals, 1, 4);
                //addback wih energy limit
                addbackv2(clover, ev, cals, 0, 0);
                addbackv2(clover, ev, cals, 1, 4);
                //modified ADDBACK
                addbackv3(clover, ev, cals, 0, 0);
                addbackv3(clover, ev, cals, 1, 4);
            }


            // Timespectra
            for(int i=0; i<ev.cnt-1; i++) {
                for(int j=i+1; j<ev.cnt; j++) {
                    // determine energy and corresponding detector-leaf
                    int enA=ev.energies[i];
                    int enB=ev.energies[j];
                    int detA=ev.dets[i];
                    int detB=ev.dets[j];

                    // energy_calibration
                    t_polynom calA=cals[detA];
                    t_polynom calB=cals[detB];

                    double x = 1.*enA;
                    double p = 1.*random()/RAND_MAX;
                    x = x + p;
                    x = calA.coeffs[0] + calA.coeffs[1]*x;
                    x = x * STRETCHFACTOR;
                    enA = (int)x;

                    x = 1.*enB;
                    p = 1.*random()/RAND_MAX;
                    x = x + p;
                    x = calB.coeffs[0] + calB.coeffs[1]*x;
                    x = x * STRETCHFACTOR;
                    enB = (int)x;

                    // determine index
                    int index=get_timespec_index(detA, detB);
                    /*if (detA==0) {
                        index=detB;
                    }
                    else if (detA==1) {
                        index=detB+2;
                    }
                    else if (detA==2) {
                        index=6;
                    }
                    else {
                        fprintf(stderr, "unexpected value pair for detectors: detA %d, detB %d\n", detA, detB);
                        continue;
                    }*/

                    // search for coincident pairs
                    //if((detA<=3)&&(detB<=3)) {
                    if((detA<LEAFCNT)&&(detB<LEAFCNT)) {
                        // trheshold for walk-correction WALK_THRESH
                        if (enA>WALK_THRESH && enB>WALK_THRESH) {

                            double q = 1.*random()/RAND_MAX;
                            double time_i=(double)ev.times[i]+q+walkcorrect(enA);

                            double r = 1.*random()/RAND_MAX;
                            double time_j=(double)ev.times[j]+r+walkcorrect(enB);

                            int64_t timediff_hist = round(time_i-time_j+(HISTCHANNELMAX/2));

                            if (!WALK_CORRECTION) {
                                timediff_hist=ev.times[i]-ev.times[j]+(HISTCHANNELMAX/2);
                                //printf("no correction, timediff: %lld\n", timediff_hist);
                            }

                            if(timediff_hist>=0&&timediff_hist<HISTCHANNELMAX) {
                                timehist[index][timediff_hist]++;
                                timehist[0][timediff_hist]++;
                            }
                        }
                    }
                }
            }


            if(createMatrices){
                for(int i=0; i<ev.cnt-1; i++) {
                    for(int j=i+1; j<ev.cnt; j++) {

                        int detA = ev.dets[i];
                        int detB = ev.dets[j];


                        int veto=0;
                        if (ev.veto[i]||ev.veto[j]) {
                            veto=1;
                        }

                        // set cut for time on coincidence spectra
                        if (TCUT) {

                            int enA=ev.energies[i];
                            int enB=ev.energies[j];

                            // energy_calibration
                            t_polynom calA=cals[detA];
                            t_polynom calB=cals[detB];

                            double x = 1.*enA;
                            double p = 1.*random()/RAND_MAX;
                            x = x + p;
                            x = calA.coeffs[0] + calA.coeffs[1]*x;
                            x = x * STRETCHFACTOR;
                            enA = (int)x;

                            x = 1.*enB;
                            p = 1.*random()/RAND_MAX;
                            x = x + p;
                            x = calB.coeffs[0] + calB.coeffs[1]*x;
                            x = x * STRETCHFACTOR;
                            enB = (int)x;

                            //determine timedifference and apply walk-correction

                            // trheshold for walk-correction
                            if (enA>WALK_THRESH && enB>WALK_THRESH) {

                                double q = 1.*random()/RAND_MAX;
                                double time_i=(double)ev.times[i]+q+walkcorrect(enA);

                                double r = 1.*random()/RAND_MAX;
                                double time_j=(double)ev.times[j]+r+walkcorrect(enB);

                                int64_t timediff = round(time_i-time_j+(HISTCHANNELMAX/2));
                                int64_t t_zero = HISTCHANNELMAX/2;

                                if ((timediff>=0&&timediff<HISTCHANNELMAX)) {
                                    if ((timediff>=(t_zero-TR))&&(timediff<=(t_zero+TR))) {
                                        add_coincidence(detA, ev.energies[i], cals[detA], detB, ev.energies[j], cals[detB], veto);
                                    }
                                }
                            }
                        }

                        // set no cut for time on coincidence spectra
                        else {
                            add_coincidence(detA, ev.energies[i], cals[detA], detB, ev.energies[j], cals[detB], veto);
                        }
                    }
                }
            }
            if(countingRate){
                int64_t diff = lasttime - firsttime;
                double ddiff = 12.5e-9*diff;
                if(total_time+ddiff>=(crTimes.size()+1)*crTimeDiff){
                    crTimes.push_back(total_time+ddiff);
                    if(crTimes.size()>2){
                        countingRates.push_back(crHits/(crTimes.at(crTimes.size()-1)-crTimes.at(crTimes.size()-2)));
                        countingRatesError.push_back(sqrt(crHits)/(crTimes.at(crTimes.size()-1)-crTimes.at(crTimes.size()-2)));
                    }
                    else{
                        countingRates.push_back(crHits/(crTimes.at(crTimes.size()-1)));
                        countingRatesError.push_back(sqrt(crHits)/(crTimes.at(crTimes.size()-1)));
                    }
                    crHits=0;
                }
            }
        }

        //Ausgabe der Laufzeiten des Sub-runs

        int64_t diff = lasttime - firsttime;
        double ddiff = 12.5e-9*diff;
        printf("\nLaufzeit des Sub-runs in Sekunden: %lf\n", ddiff);
        total_time=total_time+ddiff;

        if(countingRate){
            if(total_time>=(crTimes.size()+1)*crTimeDiff){
                crTimes.push_back(total_time);
                countingRates.push_back(crHits);
                crHits=0;
            }
        }


    }
    printf("Gesamtlaufzeit in Sekunden: %lf\n", total_time);

    printf("\n");
    /*for(int nr=0; nr<5; nr++) {
        gettestevent(&ev);
        for(int i=0; i<ev.cnt-1; i++) {
            for(int j=i+1; j<ev.cnt; j++) {
                add_coincidence(ev.dets[i], ev.energies[i], cals[i], ev.dets[j], ev.energies[j], cals[j]);
            }
        }
    }*/

    qDebug() << "Writing spectra";
    // write spectra
    for(int i=0; i<HISTMAX; i++) {
        // single
        char name[2000];
        sprintf(name, "%s%d", HISTPREFIX, i);
        write_hist(name, hists[i], HISTCHANNELMAX);
/*
        // bgo supressed
        sprintf(name, "%s%d", HISTPREFIX_BGO, i);
        write_hist(name, hists_bgo[i], HISTCHANNELMAX);
*/
        // triple
        sprintf(name, "%s%d", HISTPREFIX_TRIPLE, i);
        write_hist(name, hists_triple[i], HISTCHANNELMAX);

        // quad
        sprintf(name, "%s%d", HISTPREFIX_QUAD, i);
        write_hist(name, hists_quad[i], HISTCHANNELMAX);
    }

    qDebug() << "Writing timespectra";
    //write timespectra
    for(int i=0; i<TIMEHISTMAX; i++) {
        char name[2000];
        //sprintf(name, "%s%02d", HISTPREFIX_TIME, i);
        //write_hist(name, timehist[i], HISTCHANNELMAX);
    }


    if(addback){
        qDebug() << "Writing Addback spectra";
      for(int j=0; j<2; j++)
      {
        char name[2000];
        switch(j)
        {
          case 0:
            sprintf(name, "%s_A", HISTPREFIX_TOTAL1);
            write_hist(name, clover[j].total[0], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_ADD1);
            write_hist(name, clover[j].addback[0], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_DIRECT1);
            write_hist(name, clover[j].direct[0], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_SINGLES1);
            write_hist(name, clover[j].singles[0], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_TOTAL2);
            write_hist(name, clover[j].total[1], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_ADD2);
            write_hist(name, clover[j].addback[1], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_DIRECT2);
            write_hist(name, clover[j].direct[1], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_SINGLES2);
            write_hist(name, clover[j].singles[1], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_TOTAL3);
            write_hist(name, clover[j].total[2], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_ADD3);
            write_hist(name, clover[j].addback[2], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_DIRECT3);
            write_hist(name, clover[j].direct[2], HISTCHANNELMAX);
            sprintf(name, "%s_A", HISTPREFIX_SINGLES3);
            write_hist(name, clover[j].singles[2], HISTCHANNELMAX);
            break;

        case 1:
          sprintf(name, "%s_B", HISTPREFIX_TOTAL1);
          write_hist(name, clover[j].total[0], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_ADD1);
          write_hist(name, clover[j].addback[0], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_DIRECT1);
          write_hist(name, clover[j].direct[0], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_SINGLES1);
          write_hist(name, clover[j].singles[0], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_TOTAL2);
          write_hist(name, clover[j].total[1], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_ADD2);
          write_hist(name, clover[j].addback[1], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_DIRECT2);
          write_hist(name, clover[j].direct[1], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_SINGLES2);
          write_hist(name, clover[j].singles[1], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_TOTAL3);
          write_hist(name, clover[j].total[2], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_ADD3);
          write_hist(name, clover[j].addback[2], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_DIRECT3);
          write_hist(name, clover[j].direct[2], HISTCHANNELMAX);
          sprintf(name, "%s_B", HISTPREFIX_SINGLES3);
          write_hist(name, clover[j].singles[2], HISTCHANNELMAX);
          break;
        }
      }
    }


    if(createMatrices){
        qDebug() << "Creating matrices";
        printf("writing matrices\n");

        // write matrices
            for(int i=0; i<MATMAX; i++) {
                    flush_stack(i);
                    char name[2000];
                    get_matrixfilename(name, i);
                    read_matrix(name, &maTrix[0][0]);

                    // create directories
                    char dirname[2000];
                    printf("c%d\n",i);
                    sprintf(dirname, "matrices/c%d", i);
                    if(access(dirname, F_OK)) {
                            if(mkdir(dirname, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
                                    fprintf(stderr, "cannot created directory matrices");
                                    exit(1);
                            }
                    }

                    // write matrix file
                    sprintf(name, "%s/c%d.mtx", dirname, i);
                    write_matrix(name, &maTrix[0][0], MATCHANNELMAX, MATCHANNELMAX);

                    // write projections
                    // x
                    long unsigned int hist[MATCHANNELMAX];
                    for(int y=0; y<MATCHANNELMAX; y++) {
                            hist[y]=0;
                            for(int x=0; x<MATCHANNELMAX; x++) {
                                    hist[y]+=maTrix[x][y];
                            }
                    }
                    sprintf(name, "%s/c%d.prx", dirname, i);
                    save_spec(name, hist, MATCHANNELMAX);

                    // y
                    for(int x=0; x<MATCHANNELMAX; x++) {
                            hist[x]=0;
                            for(int y=0; y<MATCHANNELMAX; y++) {
                                   hist[x]+=maTrix[x][y];
                            }
                    }
                    sprintf(name, "%s/c%d.pry", dirname, i);
                    save_spec(name, hist, MATCHANNELMAX);

            // delete temporary matrix
            close(matrix_fds[i]);
            get_matrixfilename(name, i);
            if(remove(name)==-1) {
                            fprintf(stderr, "warning: cannot remove matrix file %s\n", name);
                            fprintf(stderr, "%s\n", strerror(errno));
                    }

            }
    }

    if(countingRate){
        QString content="";
        QFile output(countingRateFile);

        for(int i=0; i<crTimes.size(); i++){
            content+=QString::number(crTimes.at(i))+"\t"+QString::number(countingRates.at(i)) + "\t"+QString::number(countingRatesError.at(i)) + "\n";
        }

        if(output.open(QIODevice::WriteOnly)){
            output.write(content.toAscii());
            output.close();
        }
    }


    delete er;

    // find biggest coinc
    int coinc_max = 0;
    for(int i=DETMAX-1; i>0&&stats_coinc[i]==0; i--) {
        coinc_max = i;
    }

    printf("coincidences (zeros included)\n");
    for(int i=1; i<coinc_max; i++) {
        printf("%d %d\n", i, stats_coinc[i]);
    }



    if(write_summary) {
        printf("detA detB timeindex matrixindex\n");
        for(int detA=0; detA<LEAFCNT; detA++) {
            for(int detB=detA+1; detB<LEAFCNT; detB++) {
                printf("%4d %5d %8d %11d\n", detA, detB, get_timespec_index(detA, detB), get_matrix_index(detA, detB));
            }
        }
    }

    return 0;
}
