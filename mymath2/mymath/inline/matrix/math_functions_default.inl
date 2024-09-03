#pragma once

#include "../../headers/matrix.h"

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS

	template<class T, size_t n, size_t m, size_t l>
	matrix<T, n, l>* matrix_multiply_standart(
		const matrix<T, n, m>& a,
		const matrix<T, m, l>& b,
		matrix<T, n, l>* c_ptr = nullptr)
	{
		if (c_ptr == nullptr)
			c_ptr = new matrix<T, n, l>;

		if (c_ptr == nullptr)
			return nullptr;

		matrix<T, n, l>::fill(*c_ptr, 0);

		matrix<T, n, l>& c = c_ptr;

		for (int i = 0; i != n; ++i)
			for (int j = 0; j != n; ++j) {
				auto tmp = a[i][j];
				for (int k = 0; k != n; ++k)
					c[i][k] += tmp * b[j][k];
			}

		return c_ptr;
	}

	template<class T, size_t n, size_t m, size_t l>
	matrix<T, n, l>* multiply(
		const matrix<T, n, m>& a,
		const matrix<T, m, l>& b,
		matrix<T, n, l>* c_ptr = nullptr)
	{
		return matrix_multiply_standart(a, b, c_ptr);
	}
	
	template<class T, size_t n>
	matrix<T, n, n>* inv(const matrix<T, n, n>& a, matrix<T, n, n>* c_ptr = nullptr);



}