#pragma once

#include "../../headers/vector.h"
#include "../../details/__expr.inl"

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS
	template<class T, size_t n>
	double norm2(const vector<T, n>& v) {
		double tmp = 0;
		for (size_t i = 0; i != n; ++i) {
			tmp += v[i] * v[i];
		}
		return tmp;
	}
	
	template<class T, size_t n>
	double norm(const vector<T, n>& v) {
		return std::sqrt(norm2(v));
	}

}