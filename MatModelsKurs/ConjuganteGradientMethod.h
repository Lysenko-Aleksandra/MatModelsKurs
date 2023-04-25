#ifndef CONJUGATE_GRADIENT_METHOD_H 
#define CONJUGATE_GRADIENT_METHOD_H 

#include "PackedMatrix.h"
#include "Vector.h"
#include "LinearSystemSolution.h"
#include "FillingMatrixFunctions.h"

Vector conjugateGradientMethod(PackedMatrix A, Vector f, Vector x0, int M, int N);

#endif
