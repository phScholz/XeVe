#ifndef POLYNOM_H
#define POLYNOM_H

#define CAL_SHIFT 0
#define CAL_POL 1

typedef struct {
	int deg;
	// coeff[3]*x^3 + ... + coeff[0]*x^0
	double coeffs[4];
} t_polynom;

extern t_polynom read_shiftfile(const char *name);
extern t_polynom read_polynomfile(const char *name);
extern double pol_getvalue(double x, t_polynom pol);
extern void pol_print(t_polynom pol);

#endif
