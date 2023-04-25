#ifndef SOLVE_H 
#define SOLVE_H

#include <vector>

#include "Steps&Grids.h"
#include "Vector.h"
#include "FillingMatrixFunctions.h"
#include "ConjuganteGradientMethod.h"

Vector solve(int x_n, int y_n,
	double hx, double hy, double (*k1)(double x, double y),
	double (*k2)(double x, double y), double (*g1)(double y),
	double (*g2)(double y), double (*g3)(double x),
	double (*g4)(double x), double lambda_1,
	double lambda_3,
	double (*f)(double x, double y), double x_lower,
	double x_upper, double y_lower, double y_upper) {
	int M = x_n;
	int N = x_n * y_n;

	std::vector<double> xSteps = getMainSteps(x_lower, x_upper,x_n);
	std::vector<double> xMain = getMainGrid(x_lower, xSteps);
	std::vector<double> xAux = getAuxiliaryGrid(x_lower, xMain);

	double hx = xSteps[1] - xSteps[0];

	std::vector<double> ySteps = getMainSteps(y_lower, y_upper, y_n);
	std::vector<double> yMain = getMainGrid(y_lower, ySteps);
	std::vector<double> yAux = getAuxiliaryGrid(y_lower, yMain);

	double hy = ySteps[1] - ySteps[0];

	PackedMatrix A = fillMatrixA(x_n, y_n, xMain, xAux, yMain, yAux, hx, hy,
		k1, k2, g1, g2, g3, g4, lambda_1, lambda_3, f);
	 
	Vector X0 = Vector(N);

	Vector g = fillVectorG(x_n, y_n, xMain, xAux, yMain, yAux, hx, hy,
		k1, k2, g1, g2, g3, g4, lambda_1, lambda_3, f);

	Vector Xk = conjugateGradientMethod(A, g, X0, M, N);
	return Xk;
}

#endif
