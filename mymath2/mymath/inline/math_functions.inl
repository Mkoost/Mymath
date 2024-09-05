#pragma once

#include "matrix.inl.h"
#include "vector.inl.h"
#include "../details/__expr.inl"
#include "../settings.h"

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS
#if defined(MYMATH_MATRIX_STATE) && defined(MYMATH_VECTOR_STATE)
	template<class T, size_t n, size_t m>
	vector<T, n>* matrix_and_vector_multiply_standart(
		const matrix<T, n, m>& a,
		const vector<T, m>& b,
		vector<T, n>* c_ptr = nullptr) {

		if (c_ptr == nullptr)
			c_ptr = new vector<T, n>;

		if (c_ptr == nullptr)
			return nullptr;

		vector<T, n>::fill(*c_ptr, 0);

		vector<T, n>& c = *c_ptr;

		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != m; ++j) {
				auto tmp = a[i][j];
				c[i] += tmp * b[j];
			}

		return c_ptr;
		

	};

	template<class T, size_t n, size_t m>
	vector<T, n>* multiply(
		const matrix<T, n, m>& a,
		const vector<T, m>& b,
		vector<T, n>* c_ptr = nullptr) {
		return matrix_and_vector_multiply_standart(a, b, c_ptr);
	}

	template<class T, size_t n, size_t m>
	vector<T, m>* vector_and_matrix_multiply_standart(
		const vector<T, n>& a,
		const matrix<T, n, m>& b,
		vector<T, m>* c_ptr = nullptr) {

		if (c_ptr == nullptr)
			c_ptr = new vector<T, m>;

		if (c_ptr == nullptr)
			return nullptr;

		vector<T, n>::fill(*c_ptr, 0);

		vector<T, n>& c = *c_ptr;

		for (size_t j = 0; j != m; ++j) {
			auto tmp = a[j];
			for (size_t k = 0; k != m; ++k)
				c[k] += tmp * b[j][k];
		}

		return c_ptr;
	};
	
	template<class T, size_t n, size_t m>
	vector<T, n>* multiply(
		const vector<T, n>& a,
		const matrix<T, n, m>& b,
		vector<T, m>* c_ptr = nullptr) {
		return vector_and_matrix_multiply_standart(a, b, c_ptr);
	}

#endif

}