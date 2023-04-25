#ifndef SOLVE_H 
#define SOLVE_H

#include <vector>

#include "Steps&Grids.h"
#include "PackedMatrix.h"
#include "Vector.h"

PackedMatrix fillMatrixA(int x_n, int y_n, std::vector<double>xMain,
	std::vector<double>xAux,
	std::vector<double> yMain,
	std::vector<double>yAux,
	double hx, double hy, double (*k1)(double x, double y),
	double (*k2)(double x, double y), double (*g1)(double y),
	double (*g2)(double y), double (*g3)(double x),
	double (*g4)(double x), double lambda_1,
	double lambda_3,
	double (*f)(double x, double y)) {
	int N = x_n * y_n;
	int M = x_n;
	Vector g(N);
	PackedMatrix A=PackedMatrix();
	double a_m, c_m, b_m,gValue = 0;
	for (int j = 0; j < y_n; j++) {
		for (int i = 0; i < x_n; i++) {
			int m = j * M + i;////переходим к одноиндексной записи

			if ((i == 0) && (j == 0)) {////точка с координатами 0,0
				c_m=hx/2*lambda_3
					+hx*(2*hy)*k2(xMain[i],yAux[j+1])
					+hy/2*lambda_1
					+hy/(2*hx)*k1(xAux[i+1],yMain[i]);
				gValue = hx / 2 * hy / 2 * f(xMain[i], yMain[j]
					+hx / 2 * g3(xMain[i])
					+ hy / 2 * g1(yMain[j]));
				g.setElem(gValue, m);
				A.putElement(c_m, m, true);
			}

			if ((j == 0) && (i > 0) && (i <= x_n - 1)) {////нижняя граница
				b_m = -hy / (2 * hx) * k1(xAux[i], yMain[j]);
				c_m = hx / hy * k2(xMain[i], yAux[j + 1])
					+ hx * lambda_3
					+ hy / (2 * hx) * k1(xAux[i], yMain[j])
					+ hy / (2 * hx) * k1(xAux[i + 1], yMain[j]);

				if (i == x_n - 1) {////точка x_n-1;0
					gValue = (hx * hy) / 2 * f(xMain[i], yMain[j]) + hx *
						g3(xMain[i])
						+ hy / (2 * hx) * k1(xAux[i + 1], yMain[j]) * g2(yMain[j]);
				}
				else {
					gValue=(hx * hy) / 2 * f(xMain[i], yMain[j]) + hx *
						g3(xMain[i]);
				}
				g.setElem(gValue, m);
				A.putElement(b_m, m - 1, true);
				A.putElement(c_m, m, false);
			}

			if (i == 0 && j > 0 && j <= y_n - 1) {////левая граница
				c_m = hx / (2 * hy) * k2(xMain[i], yAux[j])
					+ hx / (2 * hy) * k2(xMain[i], yAux[j + 1])
					+ hy / hx * k1(xAux[i + 1], yMain[j])
					+ hy * lambda_1;
				a_m = -hx / (2 * hy) * k2(xMain[i], yAux[j]);
				if (j == y_n - 1) {///точка 0;y_n - 1
					gValue= hx / 2 * hy * f(xMain[i], yMain[j])
						+ hy * g1(yMain[j])
						+ hx / (2 * hy) * k2(xMain[i], yAux[j + 1]) * g4(xMain[i]);
				}
				else {
					gValue= hx / 2 * hy * f(xMain[i], yMain[j]) + hy *
						g1(yMain[j]);
				}
				g.setElem(gValue, m);
				A.putElement(a_m, m - x_n, true);
				A.putElement(c_m, m, false);
			}

			if (i > 0 && i <= x_n - 1 && j > 0 && j <= y_n - 1) {////внутренние точки
				c_m = hx / hy * k2(xMain[i], yAux[j])
					+ hx / hy * k2(xMain[i], yAux[j + 1])
					+ hy / hx * k1(xAux[i], yMain[j])
					+ hy / hx * k1(xAux[i + 1], yMain[j]);
				b_m = -hy / hx * k1(xAux[i], yMain[j]);
				a_m = -hx / hy * k2(xMain[i], yAux[j]);
				if ((i == x_n - 1) && (j == y_n - 1)) {///точка x_n - 1;y_n - 1
					gValue= hx * hy * f(xMain[i], yMain[j])
						+ hy / hx * k1(xAux[i + 1], yMain[j]) * g2(yMain[j])
						+ hx / hy * k2(xMain[i], yAux[j + 1]) * g4(xMain[i]);
				}
				else if (i == x_n - 1) {///правая граница
					gValue= hx * hy * f(xMain[i], yMain[j]) 
						+ hy / hx * k1(xAux[i+ 1], yMain[j]) * g2(yMain[j]);
				}
				else if (j == y_n - 1) {//верхняя граница
					gValue=hx * hy * f(xMain[i], yMain[j]) 
						+ hx / hy * k2(xMain[i],
						yAux[j + 1]) * g4(xMain[i]);
				}
				else {
					gValue = hx * hy * f(xMain[i], yMain[j]);
				}
				g.setElem(gValue, m);
				A.putElement(a_m, m - x_n, true);
				A.putElement(b_m, m - 1, false);
				A.putElement(c_m,  m, false);
			}
		}

	}
	A.putLastRowIndex(A.getA().size());
	return A;
}
PackedMatrix fillMatrixL(PackedMatrix& A, int M, int N) {
	PackedMatrix L(A);
	double a_m;
	double c_m;
	double b_m;

	double ab;
	double c;

	double a;
	double b;

	for (int i = 1; i < L.getIR().size(); i++) {///для каждого ряда
		int start = L.getIR()[i - 1];
		int finish = L.getIR()[i];
		int nonNullableElementsAmount = finish - start;///находим число ненулевых элементов

		switch (nonNullableElementsAmount)
		{
		case 1:
			c_m = std::sqrt(L.getA()[start]);
			L.setAElem(c_m, start);

			b_m = L.getA()[L.getIR()[i]] / c_m;
			a_m = L.getA()[L.getIR()[i - 1 + M]] / c_m;

			L.setAElem(b_m, L.getIR()[i]);
			L.setAElem(a_m, L.getIR()[i - 1 + M]);

			break;
		case 2:
			ab = L.getA()[start];
			c= L.getA()[start + 1];
			c_m = std::sqrt(c - ab * ab);
			L.setAElem(c_m, start + 1);
			if (i + M < N) {
				a_m = L.getA()[L.getIR()[i - 1 + M]] / c_m;
				A.setAElem(a_m, L.getIR()[i - 1 + M]);
			}
			b_m = L.getA()[L.getIR()[i]] / c_m;
			if ((i - 1) % M == 0) {
				L.setAElem(b_m, L.getIR()[i] + 1);
			}
			else if ((i - 1) % (M - 1) != 0) {
				L.setAElem(b_m, L.getIR()[i]);
			}
			break;
		case 3:
			a = L.getA()[start];
			b = L.getA()[start + 1];
			c_m = std::sqrt(L.getA()[start + 2] - a * a - b * b);
			L.setAElem(c_m, start + 2);
			if (((i - 1) >= N - M && i < N) ||
				((i - 1) < N - M) && (i % M != 0)) {
				b_m = L.getA()[L.getIR()[i + 1]] / c_m;
				L.setAElem(b_m, L.getIR()[i] + 1);
			}
			if ((i - 1) < N - M) {
				a_m = L.getA()[L.getIR()[i - 1 + M]] / c_m;
				A.setAElem(a_m, L.getIR()[i - 1 + M]);
			}
			break;
		}
	}
	return L;
}
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
