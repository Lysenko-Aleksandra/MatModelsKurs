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
	Vector second_vector(std::vector < double>{2.0, 2.0, 2.0});

	Vector summ = first_vector + second_vector;
	std::cout << summ;
	Vector mult_2 = first_vector * 2;
	std::cout << mult_2;
	Vector subs = first_vector - second_vector;
	std::cout << subs;
	double abs = first_vector.abs();
	std::cout << abs;
}
#endif
