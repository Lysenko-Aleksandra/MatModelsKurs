#include "ConjuganteGradientMethod.h"

Vector conjugateGradientMethod(PackedMatrix A, Vector f, Vector x0, int M, int N) {
	PackedMatrix L = fillMatrixL(A, M, N);
	int iterCount = 0;

	Vector r0 = f - x0 * A;
	Vector w0 = solveLinearSystem(L, r0, N, M);

	Vector sk = w0;

	double a, b;

	while (true) {
		iterCount++;
		a = (w0 * r0) / (sk * A * sk);//+
		Vector xk = x0 + sk * a;
		Vector rk = r0 - sk * A *a;

		double eps = 0.0001;
		///std::cout << rk.abs() << '\n'<< f.abs()<<'\n';
		if (rk.abs() / f.abs() < eps) {
			std::cout << "����� �������� " << iterCount;
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