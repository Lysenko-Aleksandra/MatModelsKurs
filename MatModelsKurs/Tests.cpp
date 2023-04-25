#include "Tests.h"

#include "Solve.h"
#include "Steps&Grids.h"

std::ostream& operator <<(std::ostream& out, Vector& v) {
	std::vector<double>vals = v.getValues();
	for (int i = 0; i < vals.size(); i++) {
		out << vals[i] << ' ';
	}
	out << '\n';
	return out;
}

std::ostream& operator <<(std::ostream& out, PackedMatrix& m) {
	std::vector<int> IC_ = m.getIC();
	std::vector<int> IR_ = m.getIR();
	std::vector<double> A_ = m.getA();
	for (int rowFirstIndex = 0; rowFirstIndex < IR_.size() - 1; rowFirstIndex++) {
		std::vector<int> IC_ = m.getIC();
		int start = m.getIR()[rowFirstIndex];
		int end = IR_[rowFirstIndex + 1];
		std::vector<double> row(IR_.size() - 1);
		for (int j = 0; j < IR_.size() - 1; j++) {
			for (int k = start; k < end; k++) {
				if (IC_[k] == j) {
					row[j] = A_[k];
				}
			}
		}
		for (int i = 0; i < row.size(); i++) {
			out << row[i] << '\t';
		}
		out << '\n';
	}
	return out;
}

void test_matrix_unpacked() {
	PackedMatrix m = PackedMatrix();
	m.putElement(13, 0, true);
	m.putElement(7, 1, false);
	m.putElement(1, 3, false);

	m.putElement(14, 1, true);
	m.putElement(8, 2, false);
	m.putElement(2, 4, false);

	m.putElement(15, 2, true);
	m.putElement(3, 5, false);

	m.putElement(16, 3, true);
	m.putElement(9, 4, false);
	m.putElement(4, 6, false);

	m.putElement(17, 4, true);
	m.putElement(10, 5, false);
	m.putElement(5, 7, false);

	m.putElement(18, 5, true);
	m.putElement(6, 8, false);

	m.putElement(19, 6, true);
	m.putElement(11, 7, false);

	m.putElement(20, 7, true);
	m.putElement(12, 8, false);

	m.putElement(21, 8, true);


	m.putLastRowIndex(21);
	std::cout << m;
}

void test_vector_manipulations() {
	Vector first_vector(std::vector < double>{2.0, 2.0, 2.0});
	Vector first_vector_2(std::vector < double>{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0});
	Vector second_vector(std::vector < double>{2.0, 2.0, 2.0});

	Vector summ = first_vector + second_vector;
	std::cout << summ;
	Vector mult_2 = first_vector * 2;
	std::cout << mult_2;
	Vector subs = first_vector - second_vector;
	std::cout << subs;
	double abs = first_vector.abs();
	std::cout << abs << '\n';

	PackedMatrix m = PackedMatrix();
	m.putElement(13, 0, true);
	m.putElement(7, 1, false);
	m.putElement(1, 3, false);

	m.putElement(14, 1, true);
	m.putElement(8, 2, false);
	m.putElement(2, 4, false);

	m.putElement(15, 2, true);
	m.putElement(3, 5, false);

	m.putElement(16, 3, true);
	m.putElement(9, 4, false);
	m.putElement(4, 6, false);

	m.putElement(17, 4, true);
	m.putElement(10, 5, false);
	m.putElement(5, 7, false);

	m.putElement(18, 5, true);
	m.putElement(6, 8, false);

	m.putElement(19, 6, true);
	m.putElement(11, 7, false);

	m.putElement(20, 7, true);
	m.putElement(12, 8, false);

	m.putElement(21, 8, true);


	m.putLastRowIndex(21);

	Vector subs_3 = first_vector_2 * m;

	std::cout << subs_3;
}



double k_1_first_test(double x, double y) {
	return 3.0;
}

double k_2_first_test(double x, double y) {
	return 3.0;
}

double lambda_1 = 3.0;

double lambda_3 = 3.0;


double g_1_first_test(double y) {
	return 9.0;
}

double g_2_first_test(double y) {
	return 3.0;
}

double g_3_first_test(double x) {
	return 9.0;
}

double g_4_first_test(double x) {
	return 3.0;
}


double f_first_test(double x, double y) {
	return 0.0;
}

double u_first_test(double x, double y) {
	return 3.0;
}



double k_1_second_test(double x, double y) {
	return 3.0 * x + y + 3;
}

double k_2_second_test(double x, double y) {
	return x + 3.0 * y + 3;
}


double g_1_second_test(double y) {
	return 9.0 * y + 9;
}

double g_2_second_test(double y) {
	return 3 * y + 6;
}

double g_3_second_test(double x) {
	return 9 * x * x + 9;
}

double g_4_second_test(double x) {
	return 3.0 * x * x + 6;
}


double f_second_test(double x, double y) {
	return -36 * x - 6 * y - 27;
}

double u_second_test(double x, double y) {
	return 3 * x * x + 3 * y + 3;
}


double k_1_third_test(double x, double y) {
	return 3.0 * x * x + y + 3;
}

double k_2_third_test(double x, double y) {
	return x + 3.0 * y * y + 3;
}


double g_1_third_test(double y) {
	return 9 * std::pow(y, 3) + 9;
}

double g_2_third_test(double y) {
	return 3 * y * y + 6;
}

double g_3_third_test(double x) {
	return  9 * std::pow(x, 2) + 9;
}

double g_4_third_test(double x) {
	return 3.0 * x * x + 6;
}


double f_third_test(double x, double y) {
	return -54 * x * x - 18 * x * y - 108 * std::pow(y, 3) - 60 * y - 18;
}

double u_third_test(double x, double y) {
	return 3 * x * x + 3 * y * y + 3;
}

void first_test() {
	for (int x_n = 4; x_n <= 512; x_n *= 2) {
		int y_n = x_n;

		int N = x_n * y_n;
		int M = x_n;

		double x_lower = 0.0;
		double x_upper = 1.0;

		double y_lower = 0.0;
		double y_upper = 1.0;


		std::vector<double> xSteps = getMainSteps(x_lower, x_upper, x_n);
		std::vector<double> xMain = getMainGrid(x_lower, xSteps);


		std::vector<double> ySteps = getMainSteps(y_lower, y_upper, y_n);
		std::vector<double> yMain = getMainGrid(y_lower, ySteps);


		Vector u_res_test = solve(x_n, y_n, k_1_first_test, k_2_first_test,
			g_1_first_test, g_2_first_test, g_3_first_test, g_4_first_test, lambda_1, lambda_3, f_first_test,
			x_lower, x_upper, y_lower, y_upper);

		Vector u_res(N);

		for (int j = 0; j <= y_n - 1; j++) {
			for (int i = 0; i <= x_n - 1; i++) {
				int k = M * j + i;
				u_res.setElem(u_first_test(xMain[i], yMain[j]), k);
			}
		}

		std::vector<double> epsilon;

		for (int i = 0; i < N; i++) {
			epsilon.push_back(std::abs(u_res.getValues()[i] - u_res_test.getValues()[i]));
		}

		double maxEps = *std::max_element(epsilon.begin(), epsilon.end());

		std::cout.scientific;
		std::cout << "N_x = " << x_n << " N_y = " << y_n << '\n';
		std::cout << "Max epsilon is " << maxEps << '\n';

	}

}