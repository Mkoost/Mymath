#pragma once

#include "matrix.inl.h"
#include "vector.inl.h"
#include "../details/__expr.inl"
#include "../settings.h"
#include <cmath>

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
		const vector<T, n>&		a,
		const matrix<T, n, m>&  b,
		vector<T, m>*		    c_ptr = nullptr) {

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
		const vector<T, n>&     a,
		const matrix<T, n, m>&  b,
		vector<T, m>*			c_ptr = nullptr) {
		return vector_and_matrix_multiply_standart(a, b, c_ptr);
	}

	template<class T, size_t n, typename tmp_T=double>
	vector<T, n>* qr_solve(
		matrix<T, n, n>&	 A,
		const vector<T, n>&  b,
		vector<T, n>*		 x_ptr	    =  nullptr, 
		matrix<T, n, n>*     triag_ptr  =  nullptr, 
		tmp_T				 zero	    =  0e-15) {

		char ptrs = 0;
		
		if (x_ptr == nullptr) {
			x_ptr = new vector<T, n>;
			ptrs = 1;
		}

		if (triag_ptr != nullptr)
			triag_ptr->copy(A);
		else triag_ptr = &A;

		
		vector<T, n>& x = *x_ptr;

		x.copy(b);

		for (size_t i = 0; i != n - 1; ++i) {
			for (size_t j = i + 1; j != n; ++j) {
				tmp_T c = (*triag_ptr)[i][i];
				tmp_T s = (*triag_ptr)[j][i];

				if (std::fabs((*triag_ptr)[i][j]) <= zero) continue;

				tmp_T l = std::sqrt((c * c) + (s * s));
				

				for (size_t k = i; k != n; ++k) {
					tmp_T tmp_1 = (*triag_ptr)[i][k];
					tmp_T tmp_2 = (*triag_ptr)[j][k];
					(*triag_ptr)[i][k] = (c * tmp_1 + s * tmp_2) / l;
					(*triag_ptr)[j][k] = (- s * tmp_1 + c * tmp_2) / l ;
				}

				tmp_T tmp_1 = x[i];
				tmp_T tmp_2 = x[j];
				x[i] = (c * tmp_1 + s * tmp_2) / l;
				x[j] = (-s * tmp_1 + c * tmp_2) / l;
				
				(*triag_ptr)[j][i] = 0;

				
			}
		}

		for (size_t i = 0; i != n - 1; ++i) 
			for (size_t j = i + 1; j != n; ++j) {
				if (std::fabs((*triag_ptr)[i][j]) <= zero) {
					if(ptrs) delete x_ptr;
					throw(std::invalid_argument("The matrix is singular"));
				}
			}

		for (size_t i = 1; i != n; ++i) {
			x[n - i] /= (*triag_ptr)[n - i][n - i];
			tmp_T tmp = x[n - i];
			for (size_t j = i + 1; j != n + 1; ++j)
				x[n - j] -= (*triag_ptr)[n - j][n - i] * tmp;
		}
		
		x[0] /= (*triag_ptr)[0][0];

		return x_ptr;
	}



#endif
}