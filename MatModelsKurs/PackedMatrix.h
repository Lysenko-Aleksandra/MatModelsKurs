#ifndef PACKED_MATRIX_H 
#define PACKED_MATRIX_H 

#include<vector>
#include <iostream>

class PackedMatrix {
public:
	PackedMatrix() = default;

	PackedMatrix(const PackedMatrix& m) {
		A_ = m.A_;
		IC_ = m.IC_;
		IR_ = m.IR_;
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

	void setAElem(double elem, int position) {
		A_[position] = elem;
	}
	
	PackedMatrix getTransposed(int N, int M) {
		PackedMatrix transposed= PackedMatrix();

		for (int col = 0; col < N; col++) {
			int count = 0;
			for (int i = 0; i < IC_.size(); i++) {
				if (IC_[i] == col) {
					int offset = 0;
					if (count == 1) {
						offset = 1;
					}
					if ((count == 2) || ((col + 1) % M == 0) && (count == 1)) {
						offset = M;
					}
					transposed.putElement(A_[i], col + offset, count == 0);
					count++;
				}
			}
		}
		transposed.putLastRowIndex(transposed.A_.size());
		return transposed;
	}

private:
	std::vector < double>A_;
	std::vector<int>IC_;
	std::vector<int>IR_;
};

std::ostream& operator <<(std::ostream& out, PackedMatrix& m);

#endif
