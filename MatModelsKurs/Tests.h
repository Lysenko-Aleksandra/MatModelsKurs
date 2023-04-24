#ifndef TESTS_H 
#define TESTS_H 

#include<vector>
#include <iostream>

#include "PackedMatrix.h"
#include "Vector.h"

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
	std::cout<<m;
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
	std::cout << abs<<'\n';

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

	Vector subs_3= first_vector_2 * m;

	std::cout << subs_3;
}

double k_1_first_test(double x, double y) {
	return 3.0;
}

double k_2_first_test(double x, double y) {
	return 3.0;
}

double lambda_1() {
	return 3.0;
}

double lambda_2() {
	return 3.0;
}

double g_1_first_test(double x, double y) {
	return 9.0;
}

double g_2_first_test(double x, double y) {
	return 3.0;
}

double g_3_first_test(double x, double y) {
	return 9.0;
}

double g_4_first_test(double x, double y) {
	return 3.0;
}


double f_first_test(double x, double y) {
	return 0.0;
}


void first_test() {
	for (int x_n = 4; x_n <= 512; x_n *= 2) {
		int y_n=x_n;
		double R = 2.0;
		double h = R / x_n;

		std::vector<double>x0(x_n);///��������� ����������� �� x(������� ������)

		std::vector<double>y0(y_n);///��������� ����������� �� y(������� ������)
		
		/*std::vector<std::vector<double>> u_res_test = solve(k_test, q_test, f_non_linear_test, nu_non_linear_test, R, i, true, u0, 1.0, j);*/
		std::vector<std::vector<double>> u_res_test(x_n, std::vector<double>(y_n));
		////����������, ����� solve ������ 

		std::vector<std::vector<double>> u_res(x_n, std::vector<double>(y_n));
		for (int k = 0; k < y_n; k++) {///������������� �������
			for (int n = 0; n < x_n; n++) {///������ ������
				u_res[n][k] = 3.0;
			}
		}
		std::vector<double> epsilon;
		for (int s = 0; s < y_n; s++) {///������������� �������
			for (int n = 0; n < x_n; n++) {///������ ������
				epsilon.push_back(abs(u_res[n][s] - u_res_test[n][s]));
			}
		}
		double maxEps = *std::max_element(epsilon.begin(), epsilon.end());

		std::cout.scientific;
		std::cout << "N_x = " << x_n << " N_y = " << y_n << '\n';
		std::cout << "Max epsilon is " << maxEps << '\n';

	}
}


double k_1_second_test(double x, double y) {
	return 3.0*x+y+3;
}

double k_2_second_test(double x, double y) {
	return x+3.0*y+3;
}


double g_1_second_test(double x, double y) {
	return -18.0*x-9.0*x*x+9;
}

double g_2_second_test(double x, double y) {
	return 3*x*x+6;
}

double g_3_second_test(double x, double y) {
	return -27.0*std::pow(y,4)+9* std::pow(y, 3)-27* std::pow(y, 2)+9;
}

double g_4_second_test(double x, double y) {
	return 3.0*y*y+6;
}


double f_second_test(double x, double y) {
	return -36*x-6*y-27;
}


double k_1_third_test(double x, double y) {
	return 3.0 * x*x + y + 3;
}

double k_2_third_test(double x, double y) {
	return x + 3.0 * y*y + 3;
}


double g_1_third_test(double x, double y) {
	return -18.0 *std::pow(x,3) + 9.0 * x * x + -18*x + 9;
}

double g_2_third_test(double x, double y) {
	return 3 * x * x + 6;
}

double g_3_third_test(double x, double y) {
	return -27.0 * std::pow(y, 4) + 9 * std::pow(y, 3) - 27 * std::pow(y, 2) + 9;
}

double g_4_third_test(double x, double y) {
	return 3.0 * y * y + 6;
}


double f_third_test(double x, double y) {
	return -54*x*x-18*x*y-108*std::pow(y,3)-60*y-18;
}

#endif
