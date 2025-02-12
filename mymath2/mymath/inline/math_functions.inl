#pragma once

#include "matrix.inl.h"
#include "vector.inl.h"
#include "dynamic_matrix.inl.h"
#include "dynamic_vector.inl.h"
#include "../details/__expr.inl"
#include "../settings.h"
#include <cmath>

#if defined(MYMATH_MATRIX_STATE) && defined(MYMATH_VECTOR_STATE)
namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS

	// переделать
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
	// переделать
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
	// переделать
	template<class T, size_t n, typename tmp_T=double>
	vector<T, n>* qr_solve(
		matrix<T, n, n>&     A,
		const vector<T, n>&  b,
		vector<T, n>*        x_ptr      =  nullptr, 
		matrix<T, n, n>*     triag_ptr  =  nullptr, 
		tmp_T                zero       =  0e-15) {

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

	// переделать
	template<class T, size_t n, typename tmp_T = double>
	vector<T, n>* gauss_solve(
		matrix<T, n, n>&   A,
		vector<T, n>&      b,
		vector<T, n>*      x_ptr = nullptr,
		matrix<T, n, n>*   A_ptr = nullptr,
		tmp_T              zero = 0e-15) {

		if (x_ptr == nullptr) {
			x_ptr = new vector<T, n>;
		}

		if (A_ptr != nullptr)
			A_ptr->copy(A);
		else A_ptr = &A;

		// вектор перестановок
		vector<size_t, n> rows;
		for (size_t i = 0; i < n; ++i) rows[i] = i;

		vector<T, n>* b_ptr = &b;
		vector<T, n>& x = *x_ptr;
		size_t main_elem_ind = 0;
		tmp_T tmp1;

		for (size_t i = 0; i < n - 1; ++i) {

			main_elem_ind = i;
			// выбор главного элемента
			for (size_t ii = i; ii < n; ++ii) {
				if ((*A_ptr)[rows[ii]][i] > (*A_ptr)[rows[i]][i]) {
					main_elem_ind = ii;
				}
			}
			std::swap(rows[i], rows[main_elem_ind]);

			// правую часть поделить
			(*b_ptr)[rows[i]] /= (*A_ptr)[rows[i]][i];

			// делим всю строку на главный элемент строки
			for (size_t j = i + 1; j < n; ++j) {
				(*A_ptr)[rows[i]][j] /= (*A_ptr)[rows[i]][i];
			}

			//вычитание по столбцу из вектора b
			for (size_t ii = i + 1; ii < n; ++ii) {
				(*b_ptr)[rows[ii]] -= (*A_ptr)[rows[ii]][i] * (*b_ptr)[rows[i]];

				// вычет из строки от элемента aii 
				for (size_t jj = i + 1; jj < n; ++jj) {
					(*A_ptr)[rows[ii]][jj] -= (*A_ptr)[rows[ii]][i] * ((*A_ptr)[rows[i]][jj]);
				}
			}
		}

		// последний элемент b поделить
		(*b_ptr)[rows[n - 1]] /= (*A_ptr)[rows[n - 1]][n - 1];

		// обратный ход
		x[rows[n - 1]] = (*b_ptr)[rows[n - 1]];
		for (size_t i = n - 1; i > 0; --i) {
			x[rows[i - 1]] = (*b_ptr)[rows[i - 1]];
			for (size_t j = i; j < n; ++j) {
				x[rows[i - 1]] -= (*A_ptr)[rows[i - 1]][j] * x[rows[j]];
			}
		}

		return x_ptr;
	}


}
#endif


#if defined(MYMATH_DYNAMIC_MATRIX_STATE) && defined(MYMATH_DYNAMIC_VECTOR_STATE)
#include "./math/dynamic_algebra.inl"
#include "./math/dynamic_eq_sys.inl"
#include "./math/dynamic_de_solvers.inl"
#endif