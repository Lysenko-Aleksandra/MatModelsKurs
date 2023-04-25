#ifndef SOLVE_H 
#define SOLVE_H

#include <vector>

#include "Vector.h"

Vector solve(int x_n, int y_n,
	double (*k1)(double x, double y),
	double (*k2)(double x, double y), double (*g1)(double y),
	double (*g2)(double y), double (*g3)(double x),
	double (*g4)(double x), double lambda_1,
	double lambda_3,
	double (*f)(double x, double y), double x_lower,
	double x_upper, double y_lower, double y_upper);

#endif
