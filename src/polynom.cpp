#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "polynom.h"
#include "advanced_strings.h"

t_polynom read_shiftfile(const char *name) {

	FILE *f = fopen(name, "r");
	if(!f) {
		fprintf(stderr, "cannot open shift file %s\n", name);
		exit(1);
	}
	
	int datalines = 0;
	char sbuf[300];
	double ch_new[2], ch_old[2];
	while(fgets(sbuf, 300, f)&&(datalines<2)) {
		if(as_is_comment(sbuf)) continue;
		if(strlen(as_trim(sbuf))==0) continue;

		sscanf(sbuf, "%lf %lf", &ch_old[datalines], &ch_new[datalines]);
		datalines++;
	}
	if(datalines!=2) {
		fprintf(stderr, "could not find 2 points in shiftfile %s\n", name);
		exit(1);
	}

	t_polynom pol;
	pol.deg = 1;
	pol.coeffs[3] = 0.;
	pol.coeffs[2] = 0.;
	pol.coeffs[1] = 1.*(ch_new[1]-ch_new[0])/(ch_old[1]-ch_old[0]);
	pol.coeffs[0] = ch_new[1]-pol.coeffs[1]*ch_old[1];

	fclose(f);

	return pol;
}

t_polynom read_polynomfile(const char *name) {
	FILE *f = fopen(name, "r");
	if(!f) {
		fprintf(stderr, "cannot open polynom file %s\n", name);
		exit(1);
	}

	t_polynom pol;
	for(int i=0; i<4; i++) {
		pol.coeffs[i] = 0.;
	}

	int datalines = 0;
	char sbuf[300];
	while(fgets(sbuf, 300, f)&&(datalines<4)) {
		if(as_is_comment(sbuf)) continue;
		if(strlen(as_trim(sbuf))==0) continue;

		sscanf(sbuf, "%lf", &pol.coeffs[datalines]);
		datalines++;
	}
	pol.deg = datalines-1;

	fclose(f);

	return pol;
}

// this function should not be invoked if speed is critical, make a local inline copy instead
double pol_getvalue(double x, t_polynom pol) {
	double res = pol.coeffs[0] + pol.coeffs[1]*x + pol.coeffs[2]*x*x + pol.coeffs[3]*x*x*x;
	return res;
}

void pol_print(t_polynom pol) {
	printf("%lfx^3 + %lfx^2 + %lfx^1 + %lfx^0\n", pol.coeffs[3], pol.coeffs[2], pol.coeffs[1], pol.coeffs[0]);
}
/*
int main(int argc, char *argv[]) {
	t_polynom pol;

	pol = read_shiftfile("bla.shift");
	pol_print(pol);

	printf("new value: %lf\n", pol_getvalue(30., pol));

	pol = read_polynomfile("bla.pol");
	pol_print(pol);
	printf("new value: %lf\n", pol_getvalue(30., pol));

}*/
