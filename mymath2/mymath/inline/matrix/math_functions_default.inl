#pragma once

#include "../../headers/matrix.h"
#include "../../details/__expr.inl"

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

		matrix<T, n, l>& c = *c_ptr;

		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != m; ++j) {
				auto tmp = a[i][j];
				for (size_t k = 0; k != l; ++k)
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
#if 0
	template<class T, size_t n>
	matrix<T, n, n>* inv(const matrix<T, n, n>& a, matrix<T, n, n>* c_ptr = nullptr);
#endif

	template<class T, size_t n, size_t m>
	matrix<T, m, n>* transpose(
		const matrix<T, n, m>& a,
		matrix<T, m, n>* b_ptr = nullptr) 
	{
		if (b_ptr == nullptr)
			b_ptr = new matrix<T, m, n>;

		if (b_ptr == nullptr)
			return nullptr;
		
		matrix<T, m, n>& b = b_ptr;

		for (size_t j = 0; j != m; ++j)
			for (size_t i = 0; i != n; ++i)
				b[i][j] = a[j][i];
		
		return b_ptr;
	}


}