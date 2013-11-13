 #include <QtGui/QApplication>
#include <QString>
#include <QDebug>
#include "xeve_window.h"
#include "evaluate.h"




int main(int argc, char *argv[]){

    if(argc==1){
        QApplication a(argc, argv);
        XeVe_Window w;
        w.show();

        return a.exec();
    }

    if(argc==2){
        QString help(argv[1]);
        if(help=="-h" || help=="--help"){
            fprintf(stderr, "usage: evaluate <file1> [<file2>]\n");
            fprintf(stderr, "\tdetector configuration: -c <configfile>\n");
			fprintf(stderr, "\taddback + Coincidence: -a/--addCoinc\n");
			fprintf(stderr, "\tdebug mode: -d/--debug\n");
			fprintf(stderr, "\tcreate no matrices: -nm/--noMatrices\n");
            fprintf(stderr, "\tprint index summary: -s\n");
            exit(1);
        }
    }

    if(argc>=3){
        evaluate *iSort = new evaluate;

			iSort->debug=false;
			iSort->addCoinc=false;
			iSort->createMatrices=true;	
			iSort->addback=false;	

		for(int i=0; i<argc; i++){			
			QString para(argv[i]);
			
			
			if(para=="-a" || para=="--addCoinc"){
				qDebug() << "addCoinc mode";
	        	iSort->createMatrices=true;
				iSort->addCoinc=true;
			}
			if(para=="-nm" || para=="--noMatrices"){
				iSort->createMatrices=false;	
				qDebug() << "no matrices";		
			}
			if(para=="-d" || para=="--debug"){
	        	iSort->debug=true;
			}
		}
        iSort->countingRate=true;
        iSort->crTop=188;
        iSort->crBottom=182;
        iSort->crTimeDiff=3600;
        iSort->countingRateFile="./CountingRates";

        return iSort->sort(argc, argv);

    }
}
