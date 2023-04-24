#ifndef SOLVE_H 
#define SOLVE_H

#include <vector>

//#include "Steps&Grids.h"
//
//
//std::vector<std::vector<double>> solve() {
//	double x_lower = 0;
//	double x_upper = 1;
//	int partsAmountX = 100;
//
//	double y_lower = 0;
//	double y_upper = 1;
//	int partsAmountY = 100;
//
//	std::vector<double> x_main_steps = getMainSteps(x_lower, x_upper, partsAmountX);
//	std::vector<double> x_aux_steps = getAuxiliarySteps(x_lower, x_main_steps);
//
//	std::vector<double> x_main_grid = getMainGrid(x_lower, x_main_steps);
//
//	std::vector<double> x_aux_minus_half_grid = getAuxiliaryMinusHalfGrid(x_lower, x_main_grid);
//	std::vector<double> x_aux_plus_half_grid = getAuxiliaryPlusHalfGrid(x_lower, x_main_grid);
//
//	std::vector<double> y_main_steps = getMainSteps(y_lower, y_upper, partsAmountY);
//	std::vector<double> y_aux_steps = getAuxiliarySteps(y_lower, y_main_steps);
//
//	std::vector<double> y_main_grid = getMainGrid(y_lower, y_main_steps);
//
//	std::vector<double> y_aux_minus_half_grid = getAuxiliaryMinusHalfGrid(y_lower, y_main_grid);
//	std::vector<double> y_aux_plus_half_grid = getAuxiliaryPlusHalfGrid(y_lower, y_main_grid);
//
//
//	//double h_equal = (R - zero) / n;//шаг основной сетки, сетка равномерная
//	//std::vector<double> mainGrid(n + 1, h_equal);///массив с шагами основной сетки
//
//	//std::vector<double> rMain(n + 1);///r по основной сетке(допускается случай с неоднородной сеткой)
//	//for (int i = 0; i < n + 1; i++) {
//	//	if (i == 0) {
//	//		rMain[i] = zero;
//	//	}
//	//	else {
//	//		rMain[i] = rMain[i - 1] + mainGrid[i];
//	//	}
//	//}
//	//std::vector<double> auxiliaryGrid(n + 1);//массив с шагами дополнительной сетки
//	//for (int i = 0; i < n + 1; i++) {
//	//	if (i == 0) {
//	//		auxiliaryGrid[i] = mainGrid[i + 1] / 2;
//	//	}
//	//	else if (i == n) {
//	//		auxiliaryGrid[i] = mainGrid[i] / 2;
//	//	}
//	//	else {
//	//		auxiliaryGrid[i] = (mainGrid[i] + mainGrid[i + 1]) / 2;
//	//	}
//	//}
//	//std::vector<double> rAuxMinusHalf(n + 1);///r -1/2
//	//for (int i = 1; i < n + 1; i++) {
//	//	rAuxMinusHalf[i] = (rMain[i] + rMain[i - 1]) / 2;
//	//}
//	//std::vector<double> rAuxPlusHalf(n + 1);///r +1/2
//	//for (int i = 0; i < n; i++) {
//	//	rAuxPlusHalf[i] = (rMain[i] + rMain[i + 1]) / 2;
//	//}
//	std::vector < std::vector <double>> A(n + 1, std::vector<double>(n + 1));///размерность зависит от числа шагов
//	std::vector<double> g(n + 1);
//	std::vector<std::vector<double>> v(n + 1, std::vector<double>(amountOfSteps + 1));///вектор, содержащий решения
//
//	for (int i = 0; i < n + 1; i++) {////инициализируем результирующей матрицы
//		v[i][0] = u0[i];
//	}
//
//	for (double step = 1; step < amountOfSteps + 1; step++) {///пробегаемся по всем переменным по времени
//		double t = (step - 1) * tau;
//		for (int i = 0; i < n + 1; i++) {///вычисляем элементы матрицы
//			if (i == 0) {
//				A[i][i] = -(rAuxPlusHalf[i] * k(rAuxPlusHalf[i], t)) / mainGrid[i + 1] - auxiliaryGrid[i] * rAuxPlusHalf[i] * q(rMain[i], t) / 2;
//				A[i][i + 1] = (rAuxPlusHalf[i] * k(rAuxPlusHalf[i], t)) / mainGrid[i];
//				g[i] = auxiliaryGrid[i] * rAuxPlusHalf[i] * f(rMain[i], t) / 2;
//			}
//			else if (i == n) {
//				A[i][i] = -(rAuxMinusHalf[i] * k(rAuxMinusHalf[i], t)) / mainGrid[i] - auxiliaryGrid[i] * rMain[i] * q(rMain[i], t);
//				A[i][i - 1] = (rAuxMinusHalf[i] * k(rAuxMinusHalf[i], t)) / mainGrid[i];
//				g[i] = auxiliaryGrid[i] * rMain[i] * f(rMain[i], t) + rMain[i] * nu(rMain[i], t);
//			}
//			else {
//				A[i][i] = -(rAuxMinusHalf[i] * k(rAuxMinusHalf[i], t)) / mainGrid[i] - (rAuxPlusHalf[i] * k(rAuxPlusHalf[i], t)) / mainGrid[i] - auxiliaryGrid[i] * rMain[i] * q(rMain[i], t);
//				A[i][i - 1] = (rAuxMinusHalf[i] * k(rAuxMinusHalf[i], t)) / mainGrid[i];
//				A[i][i + 1] = (rAuxPlusHalf[i] * k(rAuxPlusHalf[i], t)) / mainGrid[i + 1];
//				g[i] = auxiliaryGrid[i] * rMain[i] * f(rMain[i], t);
//			}
//		}
//
//		std::vector<std::vector<double>> inverceD(n + 1, std::vector<double>(n + 1));
//
//		for (int i = 0; i < n + 1; i++) {///инициализируем обратную D матрицу
//			if (i == 0) {
//				inverceD[i][i] = 2 / (rAuxPlusHalf[i] * auxiliaryGrid[i]);
//			}
//			else {
//				inverceD[i][i] = 1 / (rMain[i] * auxiliaryGrid[i]);
//			}
//		}
//
//		///находим B и C
//		std::vector<std::vector<double>> B = matrixMultiplication(inverceD, A);
//		std::vector<double> C = matrixAndColumnMultiplication(inverceD, g);
//
//
//		if (isExpicit) {///если необходимо найти решение явным методом Эйлера
//			std::vector<double> vMinusOneColumn(n + 1);///инициализируем i-1 ый столбец v как вектор
//			std::vector<double> v_i;
//			for (int j = 0; j < A.size(); j++) {
//				vMinusOneColumn[j] = v[j][step - 1];
//			}
//			v_i = explicitEulerMethod(B, C, tau, vMinusOneColumn);
//
//			for (int j = 0; j < A.size(); j++) {////инициализируем i тый столбец результирующей матрицы
//				v[j][step] = v_i[j];
//			}
//		}
//		else {///если необходимо найти решение неявным методом Эйлера
//			std::vector<double> vMinusOneColumn(n + 1);///инициализируем i-1 ый столбец v как вектор
//			std::vector<double> v_i;
//			for (int j = 0; j < A.size(); j++) {
//				vMinusOneColumn[j] = v[j][step - 1];
//			}
//			v_i = implicitEulerMethod(B, C, tau, vMinusOneColumn);
//
//			for (int j = 0; j < A.size(); j++) {////инициализируем i тый столбец результирующей матрицы
//				v[j][step] = v_i[j];
//			}
//		}
//	}
//	return v;
//}
#endif
