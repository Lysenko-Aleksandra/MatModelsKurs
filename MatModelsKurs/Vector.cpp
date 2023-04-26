#include "Vector.h"

std::ostream& operator <<(std::ostream& out, Vector& v) {
	std::vector<double>vals = v.getValues();
	for (int i = 0; i < vals.size(); i++) {
		out << vals[i] << ' ';
	}
	out << '\n';
	return out;
}
