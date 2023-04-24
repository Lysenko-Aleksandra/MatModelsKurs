#ifndef PACKED_MATRIX_H 
#define PACKED_MATRIX_H 

#include<vector>
#include <iostream>

class PackedMatrix {
public:
	PackedMatrix() {
		A_ = std::vector<double>();
		IC_ = std::vector<int>();
		IR_ = std::vector<int>();
	}
	
	PackedMatrix(std::vector < double>A, std::vector<int>IC, std::vector<int>IR) {
		A_ = A;
		IC_ = IC;
		IR_ = IR;
	}

	void putElement(double element, int j, bool isNewRowBegin) {
		A_.push_back(element);
		IC_.push_back(j);

		if (isNewRowBegin) {
			IR_.push_back(A_.size() - 1);
		}
	}

	void putLastRowIndex(int i) {
		IR_.push_back(i);
	}

	double getElement(int row, int column) {
		double value = 0;

		int rowStarts = IR_[row];
		int rowEnds = IR_[row + 1];

		for (int i = rowStarts; i < rowEnds; i++) {
			if (IC_[i] == column) {
				value = A_[i];
			}
		}
		return value;
	}

	std::vector<double> getA() {
		return A_;
	}

	std::vector<int> getIR() {
		return IR_;
	}

	std::vector<int> getIC() {
		return IC_;
	}

private:
	std::vector < double>A_;
	std::vector<int>IC_;
	std::vector<int>IR_;
};

std::ostream& operator <<(std::ostream& out, PackedMatrix &m) {
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

#endif
