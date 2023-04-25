#ifndef VECTOR_H 
#define VECTOR_H  

#include<vector>
#include <cmath>
#include <iostream>
#include "PackedMatrix.h"

class Vector {
public:
	Vector() = default;

	Vector(int n) {
		values_ = std::vector<double>(n,0);
	}

	Vector(std::vector<double> vals) {
		values_ = vals;
	}
	void setElem(double elem, int position) {
		values_[position] = elem;
	}

	std::vector<double> getValues() {
		return values_;
	}

	Vector operator +(const Vector& right_vector) {
		int sizeA = values_.size();
		int sizeB = right_vector.values_.size();
		std::vector<double> result(sizeA);
		if (sizeA == sizeB)
		{
			for (int i = 0; i < sizeA; i++) {
				result[i] = values_[i] + right_vector.values_[i];
			}
		}
		return Vector(result);
	}

	Vector operator -(const Vector& right_vector) {
		int sizeA = values_.size();
		int sizeB = right_vector.values_.size();
		std::vector<double> result(sizeA);
		if (sizeA == sizeB)
		{
			for (int i = 0; i < sizeA; i++) {
				result[i] = values_[i] - right_vector.values_[i];
			}
		}
		return Vector(result);
	}


	Vector operator *(const double& value) {
		std::vector<double> result(values_.size());
		for (int i = 0; i < values_.size(); i++) {
			result[i] = values_[i] * value;
		}
		
		return Vector(result);
	}

	Vector operator *(PackedMatrix& matrix) {
		int N = matrix.getIR().size()-1;
		std::vector<double> result(N);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double other = 0;
				if (j != i) {
					other = matrix.getElement(j, i) * values_[j];
				}
				result[i] += matrix.getElement(i, j) * values_[j] + other;
			}
		}

		return Vector(result);
	}

	double operator *(const Vector& right_vector) {
		int sizeA = values_.size();
		int sizeB = right_vector.values_.size();
		double result = 0;
		if (sizeA == sizeB)
		{
			for (int i = 0; i < sizeA; i++) {
				result += values_[i] * right_vector.values_[i];
			}
		}
		return result;
	}

	double abs() {
		Vector it = *this;

		return std::sqrt(it * it);
	}

	void pushBack(double elem) {
		values_.push_back(elem);
	}

private:
	std::vector<double> values_;
};


#endif