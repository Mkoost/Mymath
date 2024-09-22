#pragma once

#include "../../headers/dynamic_vector.h"
#include "../../details/__expr.inl"

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS


	template<class T>
	double norm2(const dynamic_vector<T>& v) {
		double tmp = 0;
		for (size_t i = 0; i != v.size(); ++i) {
			tmp += v[i] * v[i];
		}
		return tmp;
	}

	template<class T>
	double norm(const dynamic_vector<T>& v) {
		return std::sqrt(norm2(v));
	}

}