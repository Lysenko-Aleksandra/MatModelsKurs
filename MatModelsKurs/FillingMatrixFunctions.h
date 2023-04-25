#ifndef FILLING_MATRIX_FUNCTIONS_H
#define FILLING_MATRIX_FUNCTIONS_H

#include "PackedMatrix.h"
#include "Vector.h"

Vector fillVectorG(int x_n, int y_n, std::vector<double>xMain,
	std::vector<double>xAux,
	std::vector<double> yMain,
	std::vector<double>yAux,
	double hx, double hy, double (*k1)(double x, double y),
	double (*k2)(double x, double y), double (*g1)(double y),
	double (*g2)(double y), double (*g3)(double x),
	double (*g4)(double x), double lambda_1,
	double lambda_3,
	double (*f)(double x, double y));

PackedMatrix fillMatrixA(int x_n, int y_n, std::vector<double>xMain,
	std::vector<double>xAux,
	std::vector<double> yMain,
	std::vector<double>yAux,
	double hx, double hy, double (*k1)(double x, double y),
	double (*k2)(double x, double y), double (*g1)(double y),
	double (*g2)(double y), double (*g3)(double x),
	double (*g4)(double x), double lambda_1,
	double lambda_3,
	double (*f)(double x, double y));

PackedMatrix fillMatrixL(PackedMatrix& A, int M, int N);
#endif
