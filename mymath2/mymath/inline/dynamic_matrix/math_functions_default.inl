#pragma once

#include "../../headers/dynamic_matrix.h"
#include "../../details/__expr.inl"

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS

	template<class T>
	dynamic_matrix<T> matrix_multiply_standart(
		const dynamic_matrix<T>& a,
		const dynamic_matrix<T>& b)
	{
		if (b.rows() != a.columns) return dynamic_matrix<T>();
		dynamic_matrix<T> c(0, a.rows(), b.columns);

		size_t n = a.rows(), m = a.columns(), l = b.columns();
		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != m; ++j) {
				auto tmp = a[i][j];
				for (size_t k = 0; k != l; ++k)
					c[i][k] += tmp * b[j][k];
			}

		return c;
	}

	template<class T>
	dynamic_matrix<T> multiply(
		const dynamic_matrix<T>& a,
		const dynamic_matrix<T>& b)
	{
		return matrix_multiply_standart(a, b);
	}
#if 0
	template<class T, size_t n>
	matrix<T, n, n>* inv(const matrix<T, n, n>& a, matrix<T, n, n>* c_ptr = nullptr);


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
#endif

}