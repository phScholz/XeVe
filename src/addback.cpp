#include <stdio.h> 
#include <stdlib.h>

#include "addback.h"


void addbackv1(ab_hist *clover, t_event ev, std::vector<t_polynom> cals, int clonr, int begin)
{
    long int erg[ev.cnt], sum=0;
    //Ermittle Multiplizitaet der Events
    int m=ab_multi(ev.cnt, ev.dets, begin);
    if(m){
        //Initialisiere Energiearray
        for(int i=0; i<ev.cnt; i++){
            erg[i]=0;
        }
        //Rekalibrierung der Energien
        for (int i=0; i<ev.cnt; i++){
            if(ev.energies[i]!=0){
                erg[i]=ab_calibrate(cals[ev.dets[i]],cals[begin],ev.energies[i]);
            }
        }
        for(int i=0; i<ev.cnt; i++){
            if((ev.dets[i] < (begin+4)) && (ev.dets[i]>=begin)){
                sum+=erg[i];
                //Inkrementieren des Directspektrums

                if(m==1){
                    clover[clonr].direct[0][erg[i]]++;
                }
                //Inkrementieren des Singlesspektrums
                clover[clonr].singles[0][erg[i]]++;

            }
        }
        //Inkrementieren des Addbackspektrums
        if(m>=2){
            clover[clonr].addback[0][sum]++;
        }
        //Inkrementieren des Totalspektrums und des Statistikarrays
        clover[clonr].total[0][sum]++;
        clover[clonr].stat_addback[0][m]++;

    }

}

void addbackv2(ab_hist *clover, t_event ev, std::vector<t_polynom> cals, int clonr, int begin){
    long int erg[ev.cnt], sum=0;
    int m=ab_multi(ev.cnt, ev.dets, begin); //Ermittle Multiplizitaet der Events
    if(m){
        //Initialisiere Energiearray
        for(int i=0; i<ev.cnt; i++){
            erg[i]=0;
        }

        //Rekalibrierung der Energien
        for (int i=0; i<ev.cnt; i++){
            if(ev.energies[i]!=0){
                erg[i]=ab_calibrate(cals[ev.dets[i]],cals[begin],ev.energies[i]);
            }
        }

        for(int i=0; i<ev.cnt; i++){
            if((ev.dets[i] < (begin+4)) && (ev.dets[i]>=begin)){
                //Es werden nur Energien aufsummiert die groesser gleich dem LIMIT sind
                if(erg[i]>=(LIMIT-cals[begin].coeffs[0])/cals[begin].coeffs[1]){
                    sum+=erg[i];
                }
                else{
                    if(erg[i]!=0){
                        clover[clonr].total[1][erg[i]]++;
                        clover[clonr].addback[1][erg[i]]++;
                    }
                }

                if(erg[i]!=0){
                    if(m==1){
                        clover[clonr].direct[1][erg[i]]++;
                    }
                    clover[clonr].singles[1][erg[i]]++;
                }

            }
        }
        if(sum!=0){
            if(m>=2){
                clover[clonr].addback[1][sum]++;
            }
            clover[clonr].total[1][sum]++;
        }
        clover[clonr].stat_addback[1][m]++;

    }

}

void addbackv3(ab_hist *clover, t_event ev, std::vector<t_polynom> cals, int clonr, int begin){
    long int *hist4=clover[clonr].singles[2];
    long int *hist3=clover[clonr].direct[2];
    long int *hist2=clover[clonr].addback[2];
    long int *hist=clover[clonr].total[2];
    int m=ab_multi(ev.cnt,ev.dets, begin);
    long int erg[m], ergi[ev.cnt];
    int num=0;
    //Initialisierung des Energie-Arrays
    for(int i=0; i<m; i++){
        erg[i]=0;
    }
    //Kalibrierung der Energie
    for (int i=0; i<ev.cnt; i++){
        ergi[i]=0;
        if(ev.energies[i]!=0){
            ergi[i]=ab_calibrate(cals[ev.dets[i]],cals[begin], ev.energies[i]);
        }
    }
    int k=0;
    for(int i=0; i<ev.cnt; i++){
        if((ev.dets[i]<(begin+4)) && (ev.dets[i]>=begin)){
            erg[k]=ergi[i];
            k++;
        }
    }
    //Inkrementierung der Direct- und Singlesspektren
    for(int i=0; i<m; i++){
        hist4[erg[i]]++;
        if(erg[i]!=0 && m==1){
            hist3[erg[i]]++;
        }
    }
    long int sum[15], max=2;
    long int *p=clover[clonr].stat_addback[2];
    switch (m){
    case 1:
        if(erg[0]!= 0) hist[erg[0]]++;
        clover[clonr].stat_addback[2][1]++;
        break;
        //---------------------------DOUBLE COINC-----------------------
    case 2:
        sum[0]=erg[0];
        sum[1]=erg[1];
        sum[2]=erg[0]+erg[1];
        ab_add2(p,hist, hist2, sum, 0, 1, 2);
        break;
        //-------------------------TRIPLE COINC------------------
    case 3:
        //fill the sum-array
        sum[0] = erg[0];
        sum[1] = erg[1];
        sum[2] = erg[2];
        sum[3] = erg[0]+erg[1];
        sum[4] = erg[0]+erg[2];
        sum[5] = erg[1]+erg[2];
        sum[6] = erg[0]+erg[1]+erg[2];
        //checking for the highest histogram bar
        max=ab_getmax(hist,sum, 7);
        //begin to play with all the cases

        for(int i=0; i<7; i++){
            if(hist[sum[i]]==hist[max]) num=i;
        }
        switch (num){
        case 0: //the histogram bar of energy 0 is the highest - so increment the histogram bar of 0
            hist[sum[0]]++;
            hist2[sum[0]]++;
            ab_add2(p, hist, hist2, sum, 1, 2, 5);
            break;

        case 1: //the histogram bar of energy 1 is the highest
            hist[sum[1]]++;
            hist2[sum[1]]++;
            ab_add2(p, hist, hist2, sum, 0, 2, 4);
            break;

        case 2: //the histogram bar of energy 2 is the highest
            hist[sum[2]]++;
            hist2[sum[2]]++;
            ab_add2(p, hist, hist2, sum, 0, 1, 3);
            break;

        case 3: //the histogram bar of energy 0 + energy 1 is the highest
            hist[sum[3]]++;
            hist[sum[2]]++;
            hist2[sum[3]]++;
            hist2[sum[2]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][2]++;
            break;

        case 4: //the histogram bar of energy 0 + energy 2 is the highest
            hist[sum[4]]++;
            hist[sum[1]]++;
            hist2[sum[4]]++;
            hist2[sum[1]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][2]++;
            break;

        case 5: //the histogram bar of energy 1 + energy 2 is the highest
            hist[sum[5]]++;
            hist[sum[0]]++;
            hist2[sum[5]]++;
            hist2[sum[0]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][2]++;
            break;

        case 6: //the histogram bar of energy 0 + energy 1 + energy 2 is the highest
            hist[sum[6]]++;
            hist2[sum[6]]++;
            clover[clonr].stat_addback[2][3]++;
            break;
        }
        break;
        //--------------------------------------QUAD COINC-------------------
    case 4:
        sum[0] = erg[0];
        sum[1] = erg[1];
        sum[2] = erg[2];
        sum[3] = erg[3];
        sum[4] = erg[0]+erg[1];
        sum[5] = erg[0]+erg[2];
        sum[6] = erg[0]+erg[3];
        sum[7] = erg[1]+erg[2];
        sum[8] = erg[1]+erg[3];
        sum[9] = erg[2]+erg[3];
        sum[10] =erg[0]+erg[1]+erg[2];
        sum[11] =erg[0]+erg[1]+erg[3];
        sum[12] =erg[0]+erg[2]+erg[3];
        sum[13] =erg[1]+erg[2]+erg[3];
        sum[14] =erg[0]+erg[1]+erg[2]+erg[3];
        //checking for the highest histogram bar
        max=ab_getmax(hist,sum, 15);

        //begin to play with all the cases
        for(int i=0; i<15; i++){
            if(hist[sum[i]]==hist[max]) num=i;
        }
        switch (num){
        case 0: //the histogram bar of energy 0 is the highest
            hist[sum[0]]++;
            sum[0]=0;
            sum[4]=0;
            sum[5]=0;
            sum[6]=0;
            sum[10]=0;
            sum[11]=0;
            sum[12]=0;
            sum[14]=0;
            ab_add3(p, hist, hist2, sum, 7, 3, 8, 2, 9, 1, 13);
            break;
        case 1: //the histogram bar of energy 1 is the highest
            hist[sum[1]]++;
            sum[1]=0;
            sum[4]=0;
            sum[7]=0;
            sum[8]=0;
            sum[10]=0;
            sum[11]=0;
            sum[13]=0;
            sum[14]=0;
            ab_add3(p, hist, hist2, sum, 5, 3, 6, 2, 9, 0, 12);
            break;
        case 2: //the histogram bar of energy 2 is the highest
            hist[sum[2]]++;
            sum[2]=0;
            sum[5]=0;
            sum[7]=0;
            sum[9]=0;
            sum[10]=0;
            sum[12]=0;
            sum[13]=0;
            sum[14]=0;
            ab_add3(p, hist, hist2, sum, 4, 3, 6, 1, 8, 0, 11);
            break;
        case 3: //the histogram bar of energy 3 is the highest
            hist[sum[3]]++;
            sum[3]=0;
            sum[6]=0;
            sum[8]=0;
            sum[9]=0;
            sum[11]=0;
            sum[12]=0;
            sum[13]=0;
            sum[14]=0;
            ab_add3(p, hist, hist2, sum, 4, 2, 5, 1, 7, 0, 10);
            break;
        case 4: //the histogram bar of energy 0 + energy 1 is the highest
            hist[sum[4]]++;
            ab_add2(p, hist, hist2, sum, 2, 3, 9);
            break;
        case 5: //the histogram bar of energy 0 + energy 2 is the highest
            hist[sum[5]]++;
            ab_add2(p, hist, hist2, sum, 1, 3, 8);
            break;
        case 6: //the histogram bar of energy 0 + energy 3 is the highest
            hist[sum[6]]++;
            ab_add2(p, hist, hist2, sum, 1, 2, 7);
            break;
        case 7: //the histogram bar of energy 1 + energy 2 is the highest
            hist[sum[7]]++;
            ab_add2(p, hist, hist2, sum, 0, 3, 6);
            break;

        case 8: //the histogram bar of energy 1 + energy 3 is the highest
            hist[sum[8]]++;
            ab_add2(p, hist, hist2, sum, 0, 2, 5);
            break;

        case 9: //the histogram bar of energy 2 + energy 3 is the highest
            hist[sum[9]]++;
            ab_add2(p, hist, hist2, sum, 0, 1, 4);
            break;

        case 10: //the histogram bar of energy 0 + energy 1 + energy 2 is the highest
            hist[sum[10]]++;
            hist[sum[3]]++;
            hist2[sum[10]]++;
            hist2[sum[3]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][3]++;
            break;

        case 11: //the histogram bar of energy 0 + energy 1 + energy 3 is the highest
            hist[sum[11]]++;
            hist[sum[2]]++;
            hist2[sum[11]]++;
            hist2[sum[2]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][3]++;
            break;

        case 12: //the histogram bar of energy 0 + energy 2 + energy 3 is the highest
            hist[sum[12]]++;
            hist[sum[1]]++;
            hist2[sum[12]]++;
            hist2[sum[1]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][3]++;
            break;

        case 13: //the histogram bar of energy 1 + energy 2 + energy 3 is the highest
            hist[sum[13]]++;
            hist[sum[0]]++;
            hist2[sum[13]]++;
            hist2[sum[0]]++;
            clover[clonr].stat_addback[2][1]++;
            clover[clonr].stat_addback[2][3]++;
            break;

        case 14: //the histogram bar of energy 0 + energy 1 + energy 2 + energy3 is the highest
            hist[sum[14]]++;
            hist2[sum[14]]++;
            clover[clonr].stat_addback[2][4]++;
            break;

        }
        break;

    }

}

void ab_hist_init(ab_hist *clover, int anzahl){
    for(int i=0; i<anzahl; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<HISTCHANNELMAX; k++){
                clover[i].total[j][k]=0;
                clover[i].addback[j][k]=0;
                clover[i].direct[j][k]=0;
                clover[i].singles[j][k]=0;
            }

            for (int k=0; k<10; k++){
                clover[i].stat_addback[j][10]=0;
            }
        }
    }
}

long int ab_calibrate(t_polynom calA, t_polynom calB, long int en){
    if(en > 0){
        double x=1.*en;
        double p=1.*random()/RAND_MAX;
        x=calA.coeffs[0] + calA.coeffs[1]*en;
        x=(x-calB.coeffs[0])/calB.coeffs[1];
        x= x + p;
        en=(int)x;
    }
    else{
        en=0;
    }
    if(en>0){
        return en;
    }
    else{
        return 0;
    }
}

long int ab_getmax(long int *hist, long int *sum, int count){
    long int max=sum[0];
    for(int i=0; i<count; i++){
        if(hist[sum[i]] >= hist[max] && sum[i]!=0){
            max=sum[i];
        }
    }
    return max;
}

void ab_add2(long int *stats, long int *hist, long int *hist2, long int *sum, int x, int y, int z){
    if(hist[sum[x]]!=0 && hist[sum[y]]!=0){
        if(hist[sum[z]]> hist[sum[x]] && hist[sum[z]] > hist[sum[y]] && sum[z]!=0){
            hist[sum[z]]++;
            hist2[sum[z]]++;
            stats[1]++;
        }
        else{
            hist[sum[x]]++;
            hist[sum[y]]++;
            hist2[sum[x]]++;
            hist2[sum[y]]++;

            stats[2]++;

        }

    }
    else{
        hist[sum[z]]++;
        hist2[sum[z]]++;
        stats[1]++;
    }
}


void ab_add3(long int *stats, long int *hist, long int *hist2, long int *sum,
                                int a, int b, int c, int d, int e, int f,int g){
    long int max = ab_getmax(hist,sum, 15);
    if(hist[sum[a]] == hist[max] && sum[a] != 0){
        hist[sum[a]]++;
        hist[sum[b]]++;
        hist2[sum[a]]++;
        hist2[sum[b]]++;
        stats[1]++;
        stats[2]++;
    }
    else{
        if(hist[sum[c]] == hist[max] && sum[c] != 0){
            hist[sum[c]]++;
            hist[sum[d]]++;
            hist2[sum[c]]++;
            hist2[sum[d]]++;
            stats[1]++;
            stats[2]++;
        }
        else{
            if(hist[sum[e]] == hist[max] && sum[e] != 0){
                hist[sum[e]]++;
                hist[sum[f]]++;
                hist2[sum[e]]++;
                hist2[sum[f]]++;
                stats[1]++;
                stats[2]++;
            }
            else{
                hist[sum[g]]++;
                hist2[sum[g]]++;
                stats[3]++;
            }
        }
    }
}


int ab_multi(int cnt,int *dets, int begin){
    int m=0;

    for(int i=0; i<cnt; i++){
        if(dets[i]>=begin){
            m++;
        }
    }

    return m;
}




/*
//---------------------------AB_LOCALMAX----------------------------
int ab_localmax(double *hist, long int i)
{	
  double m1=0,m2=0;
  
  m1=(double) 1./4.*((hist[i+4]-hist[i])/4 + (hist[i+3]-hist[i])/3 + (hist[i+2]-hist[i])/2 + (hist[i+1]-hist[i]));
  m2=(double) 1./4.*((hist[i]-hist[i-4])/4 + (hist[i]-hist[i-3])/3 + (hist[i]-hist[i-2])/2 + (hist[i]-hist[i-1]));
    
  if(m1<0 && m2>0)
  {
     return 1;
  }
  else
  {
     return 0;
  }
}

//---------------------------AB_LINEAR_FUNC----------------------------
double ab_linear_func(double y1, double y2, long int x1, long int x2, long int i)
{
  double m =(double) (y2-y1)/(x2-x1);
  double n = y1 - m*x1;
  double f = n + m*i;
  return f;
}

//---------------------------AB_MIN_LEFT----------------------------
long int ab_min_left(double *histogram, long int channel, int range)
{
  long int minimum=histogram[channel];
  long int minchannel=channel;
  
  for(int i=0; i<=range; i++)
  {
    if(histogram[channel-i]<minimum)
    {
      minchannel=channel-i;
    }
  }
  return minchannel;
}

//---------------------------AB_MIN_RIGHT----------------------------
long int ab_min_right(double *histogram, long int channel, int range)
{
  long int minimum=histogram[channel];
  long int minchannel=channel;
  
  for(int i=0; i<=range; i++)
  {
    if(histogram[channel+i]<minimum)
    {
      minchannel=channel+i;
    }
  }
  return minchannel;
}

//---------------------------AB_PEAKCHECK----------------------------
int ab_peakcheck(long int *hist, long int channel)
{
  
  //define an array for the rolling mean of the hist
  double rmean[2][AB_HISTMAX];
  //set all entries on 0  
  for(long int i=0; i<AB_HISTMAX; i++)
  {
    rmean[0][i]=0;
    rmean[1][i]=0;
  }
  //creating the rolling mean; every entry is the average of AB_WIDTH histbin's
  for(long int i=0; i<AB_HISTMAX-AB_WIDTH; i++)
  {
   for(long int j=i; j<i+AB_WIDTH; j++)
   {
      rmean[1][i]+=(double) hist[j]/AB_WIDTH;
   }  
  }
  long int a=AB_HISTMAX-AB_WIDTH/2;
  
  //shift the rolling mean in the right position
  for(long int i=-AB_WIDTH/2; i<AB_HISTMAX; i++)
  {
    rmean[0][i+AB_WIDTH]=rmean[1][i];
    rmean[1][i]=0;
  }

  
  //checking for local maximas, if ab_localmax=1 writing histvalue of bin into the array, otherwise write 0
  for(long int i=4; i<a-4;i++)
  {
    if(ab_localmax(rmean[0],i))
    {
      rmean[1][i]=hist[i];
    }
    else
    {
      rmean[1][i]=0;      
    }    
  }

  for(long int i=0; i<a;i++)
   {
      if(rmean[1][i]>0)
      {
	long int minleft=ab_min_left(rmean[0],i,AB_MINRANGE);
	long int minright=ab_min_right(rmean[0],i,AB_MINRANGE);
      
	for(long int j=i-(i-minleft); j<i+(minright-1); j++)
	{
	  if(rmean[0][j] > ab_linear_func(rmean[0][minleft],rmean[0][minright], minleft, minright, j))
	  {
	    rmean[1][j]=hist[j];
	  }
	  else {rmean[1][j]=0;}
	 }
	}
   }
  if(rmean[1][channel] > 0) return 1;
  else return 0;
}

//---------------------------AB_CALIBRATE----------------------------
long int ab_calibrate(t_polynom calA, t_polynom calB, long int en, int det)
{
  if(det!=0)
  {
    double x=1.*en;
    double p=1.*random()/RAND_MAX;
    x=calA.coeffs[0] + calA.coeffs[1]*en;
    x=(x-calB.coeffs[0])/calB.coeffs[1];
    x= x + p;
    en=(int)x;
  } 
  return en;
}
//---------------------------AB_GETMAX----------------------------
long int ab_getmax(long int *hist, long int *sum, int count)
{
  long int max=sum[0];
  for(int i=0; i<count; i++)
  {
    if(hist[sum[i]] >= hist[max] && sum[i]!=0)
    {
      max=sum[i];
    }
  }
  return max;
}

void ab_add2(long int *hist, long int *sum, long int x, long int y, long int z)
{
      if(hist[sum[z]]> hist[sum[x]] && hist[sum[z]] > hist[sum[y]])
      {
	hist[sum[z]]++;
      }
      else
      {
	hist[sum[x]]++;
	hist[sum[y]]++;	
      }
}

//---------------------------ADDBACK----------------------------
void addback(long int *hist, t_event ev, t_polynom *cals)
{
  long int erg[ev.cnt]; 
  int num=0, hits=ev.cnt;
  
  for (int i=0; i<ev.cnt; i++)
  {
    if(ev.energies[i]!=0)
    {
      erg[i]=ab_calibrate(cals[0],cals[i], ev.energies[i], ev.dets[i]);    
    }
    else
    {
      hits-=1;
    }
  }
  long int sum[15],max=0;
  
  switch (hits)
  {
    case 1:
      if(erg[0] != 0) hist[erg[0]]++; break;
      
    //---------------------------DOUBLE COINC-----------------------
    case 2:
      for(int i=0; i<hits; i++)
      {
	sum[i]=erg[i];
      }
      sum[ev.cnt]=erg[0]+erg[1];
      
      for(int i=0; i<2; i++)
      {
	if(!ab_peakcheck(hist, sum[i]))
	{
	  sum[i]=0;
	}
      }
      max=ab_getmax(hist,sum,3);
      if(hist[sum[2]]==hist[max])
      {
	hist[sum[2]]++;
      }
      else
      {
	hist[sum[0]]++;
	hist[sum[0]]++;	
      }
      break;
    //-------------------------TRIPLE COINC------------------  
    case 3:
      //fill the sum-array
      sum[0] = erg[0];
      sum[1] = erg[1];
      sum[2] = erg[2];
      sum[3] = erg[0]+erg[1];
      sum[4] = erg[0]+erg[2];
      sum[5] = erg[1]+erg[2];
      sum[6] = erg[0]+erg[1]+erg[2];
      
      //checking if there is a peak or not - if not, the entry's in the sum-array will be deleted
      for(int i=0; i<6; i++)
      {
	if(!ab_peakcheck(hist, sum[i])) 
	{
	  sum[i]=0;
	}
      }
      
      //checking for the highest histbin
      max=ab_getmax(hist,sum, 7);
      
      //begin to play with all the cases
      for(int i=0; i<7; i++)
      {
	if(hist[sum[i]]==max) num=i;
      }
      
      switch (num)
      {
	case 0: //the histbin of energy 0 is the highest 
	  hist[sum[0]]++;
	  sum[0]=0;
	  sum[3]=0;
	  sum[4]=0;
	  sum[6]=0;
	  max = ab_getmax(hist,sum,7);
      ab_add2(hist, sum, 1, 2, 5, max);
	  break;
	  
	case 1: //the histbin of energy 1 is the highest 
	  hist[sum[1]]++; 
	  sum[1]=0;
	  sum[3]=0;
	  sum[5]=0;
	  sum[6]=0;
	  max = ab_getmax(hist,sum,7);
      ab_add2(hist, sum, 0, 2, 4, max);
	  break;
	  
	case 2: //the histbin of energy 2 is the highest 
	  hist[sum[2]]++; 
          sum[2]=0;
	  sum[4]=0;
	  sum[5]=0;
	  sum[6]=0;
	  max = ab_getmax(hist,sum,7);
      ab_add2(hist, sum, 0, 1, 3, max);
	  break;
	  
	case 3: //the histbin of energy 0 + energy 1 is the highest 
	  hist[sum[3]]++;
	  hist[sum[2]]++;
	  break;
	  
	case 4: //the histbin of energy 0 + energy 2 is the highest 
	  hist[sum[4]]++;
	  hist[sum[1]]++;
	  break;
	case 5: //the histbin of energy 1 + energy 2 is the highest 
	  hist[sum[5]]++;
	  hist[sum[0]]++;
	  break;
	case 6: //the histbin of energy 0 + energy 1 + energy 2 is the highest 
	  hist[sum[6]]++;
	  break;
      }
            
      break;
      
    //--------------------------------------QUAD COINC------------------------------------
    case 4:
	sum[0] = erg[0];
	sum[1] = erg[1];
	sum[2] = erg[2];
	sum[3] = erg[3];
	sum[4] = erg[0]+erg[1];
	sum[5] = erg[0]+erg[2];
	sum[6] = erg[0]+erg[3];
	sum[7] = erg[1]+erg[2];
	sum[8] = erg[1]+erg[3];
	sum[9] = erg[2]+erg[3];
	sum[10] = erg[0]+erg[1]+erg[2];
	sum[11] = erg[0]+erg[1]+erg[3];
	sum[12] = erg[0]+erg[2]+erg[3];
	sum[13] = erg[1]+erg[2]+erg[3];
	sum[14] = erg[0]+erg[1]+erg[2]+erg[3];

	//checking if there is a peak or not - if not, the entry in the sum-array will be deleted
	for(int i=0; i<14; i++)
	{
	  if(!ab_peakcheck(hist, sum[i]))
	  {
	    sum[i]=0;
	  }
	}
	
	//checking for the highest histbin
	max=ab_getmax(hist,sum, 15); 
	
	//begin to play with all the cases
	for(int i=0; i<15; i++)
	{
	  if(hist[sum[i]]==max) num=i;
	}
	
	switch (num)
	{
	  case 0: //the histbin of energy 0 is the highest 
	    hist[sum[0]]++;
	    sum[0]=0;
	    sum[4]=0;
	    sum[5]=0;
	    sum[6]=0;
	    sum[10]=0;
	    sum[11]=0;
	    sum[12]=0;
	    sum[14]=0;
	    
	    max = ab_getmax(hist,sum, 15);
	    
	    if(hist[sum[7]] == max)
	    {
	      hist[sum[7]]++;
	      hist[sum[3]]++;
	    }
	    else
	    {
	      if(hist[sum[8]] == max)
	      {	
		hist[sum[8]]++;
		hist[sum[2]]++;
	      }
	      else
	      {
		if(hist[sum[9]] == max)
		{
		  hist[sum[9]]++;
		  hist[sum[1]]++;
		}
		else
		{
		  if(hist[sum[13]] == max)
		  {
		    hist[sum[13]]++;		    
		  }
		}
	      }
	    }
	    break;
	    
	  case 1: //the histbin of energy 1 is the highest 
	    hist[sum[1]]++;
	    sum[1]=0;
	    sum[4]=0;
	    sum[7]=0;
	    sum[8]=0;
	    sum[10]=0;
	    sum[11]=0;
	    sum[13]=0;
	    sum[14]=0;
	    
	    max = ab_getmax(hist,sum, 15);
	    if(hist[sum[5]] == max)
	    {
	      hist[sum[5]]++;
	      hist[sum[3]]++;
	    }
	    else
	    {
	      if(hist[sum[6]] == max)
	      {	
		hist[sum[6]]++;
		hist[sum[2]]++;
	      }
	      else
	      {
		if(hist[sum[9]] == max)
		{
		  hist[sum[9]]++;
		  hist[sum[0]]++;
		}
		else
		{
		  if(hist[sum[12]] == max)
		  {
		    hist[sum[12]]++;		    
		  }
		}
	      }
	    }
	    break;
	    
	  case 2: //the histbin of energy 2 is the highest 
	    hist[sum[2]]++;
	    sum[2]=0;
	    sum[5]=0;
	    sum[7]=0;
	    sum[9]=0;
	    sum[10]=0;
	    sum[12]=0;
	    sum[13]=0;
	    sum[14]=0;
  
	    max = ab_getmax(hist,sum, 15);
	    if(hist[sum[4]] == max)
	    {
		hist[sum[4]]++;
		hist[sum[3]]++;
	    }
	    else
	    {
	      if(hist[sum[6]] == max)
	      {
		hist[sum[6]]++;
		hist[sum[1]]++;
	      }
	      else
	      {
		if(hist[sum[8]] == max)
		{
		  hist[sum[8]]++;
		  hist[sum[0]]++;
		}
		else
		{
		  if(hist[sum[11]] == max)
		  {
		    hist[sum[11]]++;
		  }
		}
	      }
	    }
	    break;
	    
	  case 3: //the histbin of energy 3 is the highest 
	      hist[sum[3]]++;
	      sum[3]=0;
	      sum[6]=0;
	      sum[8]=0;
	      sum[9]=0;
	      sum[11]=0;
	      sum[12]=0;
	      sum[13]=0;
	      sum[14]=0;

	      max = ab_getmax(hist,sum, 15);
	      if(hist[sum[4]] == max)
	      {
		hist[sum[4]]++;
		hist[sum[2]]++;
	      }
	      else
	      {
		if(hist[sum[5]] == max)
		{
		  hist[sum[5]]++;
		  hist[sum[1]]++;
		}
		else
		{
		  if(hist[sum[7]] == max)
		  {
		    hist[sum[7]]++;
		    hist[sum[0]]++;
		  }
		  else
		  {
		    if(hist[sum[10]] == max)
		    {
		      hist[sum[10]]++;
		    }
		   }
		 }
	      }
	    break;

	  case 4: //the histbin of energy 0 + energy 1 is the highest 
	    hist[sum[4]]++; 
	    sum[0]=0;
	    sum[1]=0;
	    sum[4]=0;
	    sum[5]=0;
	    sum[6]=0;
	    sum[7]=0;
	    sum[8]=0;
	    sum[10]=0;
	    sum[11]=0;
	    sum[12]=0;
	    sum[13]=0;
	    sum[14]=0;
	    max = ab_getmax(hist,sum,15);
        ab_add2(hist, sum, 2, 3, 9, max);

	  case 5: //the histbin of energy 0 + energy 2 is the highest 
	    hist[sum[5]]++; 
	    if(hist[sum[8]] > hist[sum[1]] && hist[sum[8]]> hist[sum[3]])
	    {
	      hist[sum[8]]++;
	    }
	    else
	    {
	      hist[sum[1]]++;
	      hist[sum[3]]++;
	    }
	    break;

	  case 6: //the histbin of energy 0 + energy 3 is the highest 
	    hist[sum[6]]++; 
	    if(hist[sum[7]] > hist[sum[1]] && hist[sum[7]]> hist[sum[2]])
	    {
	      hist[sum[7]]++;
	    }
	    else
	    {
	      hist[sum[1]]++;
	      hist[sum[2]]++;
	    }
	    break;

	  case 7: //the histbin of energy 1 + energy 2 is the highest 
	    hist[sum[7]]++; 
	    if(hist[sum[6]] > hist[sum[0]] && hist[sum[6]]> hist[sum[3]])
	    {
	      hist[sum[6]]++;
	    }
	    else
	    {
	      hist[sum[0]]++;
	      hist[sum[3]]++;
	    }
	    break;

	  case 8: //the histbin of energy 1 + energy 3 is the highest 
	    hist[sum[8]]++; 
	    if(hist[sum[5]] > hist[sum[0]] && hist[sum[5]]> hist[sum[2]])
	    {
	      hist[sum[5]]++;
	    }
	    else
	    {
	      hist[sum[0]]++;
	      hist[sum[2]]++;
	    }
	    break;

	  case 9: //the histbin of energy 2 + energy 3 is the highest 
	    hist[sum[9]]++; 
	    if(hist[sum[4]] > hist[sum[1]] && hist[sum[4]]> hist[sum[0]])
	    {
	      hist[sum[4]]++;
	    }
	    else
	    {
	      hist[sum[0]]++;
	      hist[sum[1]]++;
	    }
	    break;
	    
	  case 10: //the histbin of energy 0 + energy 1 + energy 2 is the highest 
	    hist[sum[10]]++;
	    hist[sum[3]]++;
	    break;
	    
	  case 11: //the histbin of energy 0 + energy 1 + energy 3 is the highest 
	    hist[sum[11]]++;
	    hist[sum[2]]++;
	    break;
	    
	  case 12: //the histbin of energy 0 + energy 2 + energy 3 is the highest 
	    hist[sum[12]]++;
	    hist[sum[1]]++;
	    break;
	    
	  case 13: //the histbin of energy 1 + energy 2 + energy 3 is the highest 
	    hist[sum[13]]++;
	    hist[sum[0]]++;
	    break;
	    
	  case 14: //the histbin of energy 0 + energy 1 + energy 2 + energy 3 is the highest 
	    hist[sum[14]]++;
	    break;
	}
    }
}

//---------------------------AB_POW----------------------------
int ab_pow(int x, int y)
{
  //x^y
  for(int i=0; i<y; i++)
  {
    x*=x;
  }
  
  return x;
}
*/
