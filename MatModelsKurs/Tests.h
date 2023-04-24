#ifndef TESTS_H 
#define TESTS_H 

#include<vector>
#include <iostream>

#include "PackedMatrix.h"

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
	m.printUnpackedMatrix(std::cout);
}

#endif
