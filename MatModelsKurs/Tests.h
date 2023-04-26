#ifndef TESTS_H 
#define TESTS_H 

#include<vector>
#include <iostream>

#include "PackedMatrix.h"
#include "Vector.h"


std::ostream& operator <<(std::ostream& out, Vector& v);

std::ostream& operator <<(std::ostream& out, PackedMatrix& m);

void test_matrix_unpacked();

void test_vector_manipulations();

void first_test();
void second_test();
void third_test();

#endif
