#ifndef SOLVE_H 
#define SOLVE_H

#include <vector>

#include "Steps&Grids.h"
#include "PackedMatrix.h"
#include "Vector.h"
#include "LinearSystemSolution.h"

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
Vector conjugateGradientMethod(PackedMatrix A, Vector f, Vector x0,int M,int N) {
	PackedMatrix L = fillMatrixL(A, M, N);
	int iterCount = 0;

	Vector r0 = f - x0* A;
	Vector w0 = solveLinearSystem(L, r0,N,M);

	Vector sk = w0;

	double a, b;

	while (true) {
		iterCount++;
		a = (w0 * r0) / (sk * A * sk);
		Vector xk = x0 + sk * a;
		Vector rk = r0 = sk * A * a;

		double eps = 0.00001;
		if (rk.abs() / f.abs() < eps) {
			std::cout << "число итераций " << iterCount;
			return xk;
		}

		Vector wk = solveLinearSystem(L, rk, N, M);
		b = (wk * rk) / (w0 * r0);
		sk = wk + sk * b;
		r0 = rk;
		w0 = wk;
		x0 = xk;
	}
}

#endif
