#include "LinearSystemSolution.h"

Vector Gauss(PackedMatrix L, Vector r, int M, int N) {
	Vector y = Vector();
	y.pushBack(r.getValues()[0] / L.getA()[0]);

	for (int i = 1; i < M; i++) {
		int start = L.getIR()[i];
		y.pushBack(r.getValues()[i] - L.getA()[start] *
			y.getValues()[i - 1] / L.getA()[start + 1]);
	}

	for (int i = M; i < N; i++) {
		int start = L.getIR()[i];

		if (i % M == 0) {
			y.pushBack(r.getValues()[i] - L.getA()[start] *
				y.getValues()[i - M] / L.getA()[start + 1]);
		}
		else {
			y.pushBack(r.getValues()[i] - L.getA()[start] *
				y.getValues()[i - M]
				- L.getA()[start + 1] * y.getValues()[i - 1] / L.getA()[start + 2]);
		}
	}
	return y;
}

Vector GaussT(PackedMatrix LT, Vector y, int M, int N) {
	Vector w = Vector(N);

	int start = LT.getIR()[N - 1];
	w.setElem(y.getValues()[N - 1] / LT.getA()[start], N - 1);
	double elem = 0;

	for (int i = N - 2; i >= (N - M); i--) {
		start = LT.getIR()[i];
		elem = y.getValues()[i]
			- LT.getA()[start + 1] * w.getValues()[i + 1] / LT.getA()[start];
		w.setElem(elem, i);
	}

	for (int i = N - M-1; i >= 0; i--) {
		start = LT.getIR()[i];

		if ((i + 1) % M == 0) {
			elem = y.getValues()[i]
				- LT.getA()[start + 1] * w.getValues()[i + M] / LT.getA()[start];
		}
		else {
			elem = y.getValues()[i]
				- LT.getA()[start + 2] * w.getValues()[i + M]
				- LT.getA()[start + 1] * w.getValues()[i + 1] / LT.getA()[start];
		}
		w.setElem(elem, i);
	}
	return w;
}

Vector solveLinearSystem(PackedMatrix L, Vector r, int N, int M) {
	Vector y = Gauss(L, r, M, N);
	PackedMatrix LT = L.getTransposed(N, M);
	Vector w= GaussT(LT, y, M, N);
	return w;
}