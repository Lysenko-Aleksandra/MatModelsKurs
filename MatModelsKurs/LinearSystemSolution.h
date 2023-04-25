#ifndef LINEAR_SYSTEM_SOLUTION_H 
#define LINEAR_SYSTEM_SOLUTION_H

#include "PackedMatrix.h"
#include "Vector.h"

Vector Gauss(PackedMatrix L, Vector r, int M, int N);

Vector GaussT(PackedMatrix LT, Vector y, int M, int N);

Vector solveLinearSystem(PackedMatrix L, Vector r, int N, int M);
#endif
