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

	template<class T>
	double oct_norm(const dynamic_vector<T>& a) {
		double mx = 0;
		for (size_t j = 0; j != a.columns(); ++j) {
			double tmp = 0;
			for (size_t i = 0; i != a.rows(); ++i)
				tmp += std::fabs(a[i][j]);
			mx = max(mx, tmp);
		}
		return mx;
	}

	template<class T>
	double sphere_norm(const dynamic_vector<T>& a) {
		double mx = 0;
		for (size_t j = 0; j != a.columns(); ++j) {
			double tmp = 0;
			for (size_t i = 0; i != a.rows(); ++i)
				tmp += a[i][j] * a[i][j];
			mx = max(mx, tmp);
		}
		return std::sqrt(mx);
	}

}