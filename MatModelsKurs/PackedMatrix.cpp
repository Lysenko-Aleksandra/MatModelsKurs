#include "PackedMatrix.h"

std::ostream& operator <<(std::ostream& out, PackedMatrix& m) {
	std::vector<int> IC_ = m.getIC();
	std::vector<double> A_ = m.getA();
	for (int rowFirstIndex = 0; rowFirstIndex < m.getIR().size() - 1; rowFirstIndex++) {
		int start = m.getIR()[rowFirstIndex];
		int end = m.getIR()[rowFirstIndex + 1];
		std::vector<double> row(m.getIR().size() - 1);
		for (int j = 0; j < m.getIR().size() - 1; j++) {
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