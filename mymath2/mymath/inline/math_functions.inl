#pragma once

#include "matrix.inl.h"
#include "vector.inl.h"
#include "dynamic_matrix.inl.h"
#include "dynamic_vector.inl.h"
#include "../details/__expr.inl"
#include "../settings.h"
#include <cmath>

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS
#if defined(MYMATH_MATRIX_STATE) && defined(MYMATH_VECTOR_STATE)
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

#endif

#if defined(MYMATH_DYNAMIC_MATRIX_STATE) && defined(MYMATH_DYNAMIC_VECTOR_STATE)

	template<class T>
	dynamic_vector<T> matrix_and_vector_multiply_standart(
		const dynamic_matrix<T>& a,
		const dynamic_vector<T>& b) {

		if (a.rows() != b.size()) throw(std::invalid_argument("Matrix and vector have different sizes"));
	
		dynamic_vector<T> c(0, b.size());
		size_t n = a.rows(), m = b.size();

		for (size_t i = 0; i != n; ++i)
			for (size_t j = 0; j != m; ++j) {
				auto tmp = a[i][j];
				c[i] += tmp * b[j];
			}

		return c;


	};

	template<class T>
	dynamic_vector<T> multiply(
		const dynamic_matrix<T>& a,
		const dynamic_vector<T>& b) {
		return matrix_and_vector_multiply_standart(a, b);
	}

	template<class T>
	dynamic_vector<T> matrix_and_vector_multiply_standart(
		const dynamic_vector<T>& a,
		const dynamic_matrix<T>& b) {

		if (b.columns() != a.size()) throw(std::invalid_argument("Matrix and vector have different sizes"));

		dynamic_vector<T> c(0, a.size());
		size_t n = b.rows(), m = a.size();

		for (size_t j = 0; j != m; ++j) {
			auto tmp = a[j];
			for (size_t k = 0; k != m; ++k)
				c[k] += tmp * b[j][k];
		}

		return c;
	};

	template<class T>
	dynamic_vector<T> multiply(
		const dynamic_vector<T>& a,
		const dynamic_matrix<T>& b) {
		return vector_and_matrix_multiply_standart(a, b);
	}

	template<class T, typename tmp_T = double>
	data_structs::base_data_dynamic_vector_matrix<T, 1, 1> qr_solve(
		const dynamic_matrix<T>&    A,
		const dynamic_vector<T>&    b,
		tmp_T                       zero = 0e-15) {


		if (A.rows() != A.columns()) throw(std::invalid_argument("Sizes of matrix rows and columns must be equal"));
		if (A.rows() != b.size()) throw(std::invalid_argument("Matrix and vector have different sizes"));

		dynamic_vector<T> x(b);


		dynamic_matrix<T> triag(A);
		size_t n = b.size();

		for (size_t i = 0; i != n - 1; ++i) {
			for (size_t j = i + 1; j != n; ++j) {
				if (std::fabs(triag[i][j]) <= zero) continue;

				tmp_T c = triag[i][i];
				tmp_T s = triag[j][i];
				tmp_T l = std::sqrt((c * c) + (s * s));


				for (size_t k = i; k != n; ++k) {
					tmp_T tmp_1 = triag[i][k];
					tmp_T tmp_2 = triag[j][k];
					triag[i][k] = (c * tmp_1 + s * tmp_2) / l;
					triag[j][k] = (-s * tmp_1 + c * tmp_2) / l;
				}

				tmp_T tmp_1 = x[i];
				tmp_T tmp_2 = x[j];
				x[i] = (c * tmp_1 + s * tmp_2) / l;
				x[j] = (-s * tmp_1 + c * tmp_2) / l;

				triag[j][i] = 0;


			}
		}

		for (size_t i = 0; i != n; ++i)
			if (std::fabs(triag[i][i]) <= zero) {
				throw(std::invalid_argument("The matrix is singular"));
			}

		for (size_t i = 1; i != n; ++i) {
			x[n - i] /= triag[n - i][n - i];
			tmp_T tmp = x[n - i];
			for (size_t j = i + 1; j != n + 1; ++j)
				x[n - j] -= triag[n - j][n - i] * tmp;
		}

		x[0] /= triag[0][0];

		data_structs::base_data_dynamic_vector_matrix<T, 1, 1> res;
		res.mat[0].move(triag);
		res.vec[0].move(x);

		
		return res;
	}


	// Ќ» ќћ” Ќ≈ ѕќ ј«џ¬ј“№ !!!!!!!!!!!!!!!!! »Ќ—”Ћ№“ √ј–јЌ“»–ќ¬јЌ !!!!!!!!!!!!!
	template<class T, typename tmp_T = double>
	data_structs::base_data_dynamic_vector_matrix<T, 1, 1> gauss_solve(
		dynamic_matrix<T>& A,
		dynamic_vector<T>& b,
		tmp_T         zero = 0e-15) {

		if (A.rows() != A.columns()) throw(std::invalid_argument("Sizes of matrix rows and columns must be equal"));
		if (A.rows() != b.size()) throw(std::invalid_argument("Matrix and vector have different sizes"));

		size_t n = b.size();

		dynamic_vector<T> x(0, n);
		dynamic_matrix<T> triag(A);
		dynamic_matrix<T>* A_ptr = &triag;

		// вектор перестановок
		dynamic_vector<size_t> rows;

		rows.move(new size_t[n], n);

		for (size_t i = 0; i < n; ++i) rows[i] = i;



		dynamic_vector<T> b_(b);
		dynamic_vector<T>* b_ptr = &b_;

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

		for (size_t i = 0; i != n; ++i)
			if (std::fabs((*A_ptr)[rows[i]][i]) <= zero) {
				throw(std::invalid_argument("The matrix is singular"));
			}

		// обратный ход
		x[rows[n - 1]] = (*b_ptr)[rows[n - 1]];
		for (size_t i = n - 1; i > 0; --i) {
			x[rows[i - 1]] = (*b_ptr)[rows[i - 1]];
			for (size_t j = i; j < n; ++j) {
				x[rows[i - 1]] -= (*A_ptr)[rows[i - 1]][j] * x[rows[j]];
			}
		}


		dynamic_vector<T> x_res(0, 4);
		for (size_t i = 0; i != n; ++i) x_res[i] = x[rows[i]];


		data_structs::base_data_dynamic_vector_matrix<T, 1, 1> res;
		res.mat[0].move((*A_ptr));
		res.vec[0].move(x_res);

		return res;
	}


	template<class T, typename tmp_T = double>
	dynamic_vector<T>
	jacobi_iteration(
		const dynamic_matrix<T>&       A,
		const dynamic_vector<T>&       b,
		tmp_T                          eps, 
		dynamic_vector<T>              start = dynamic_vector<T>()){

		if (A.rows() != A.columns()) throw(std::invalid_argument("Sizes of matrix rows and columns must be equal"));
		if (A.rows() != b.size()) throw(std::invalid_argument("Matrix and vector have different sizes"));

		if (start.size() == 0) start.copy(b);
		else if (start.size() != b.size()) throw(std::invalid_argument("Vector and vector have different sizes"));

		size_t n = b.size();

		dynamic_vector<T> x(b);
		tmp_T nrm = 0;
		size_t lp = 0;
		do{
			++lp;
			nrm = 0;

			for (size_t j = 0; j != n; ++j)
				for (size_t i = 0; i != n; ++i) {
					if (i == j) continue;
					x[j] -= A[j][i] * start[i];
				}
		
			for (size_t i = 0; i != n; ++i){
				nrm = std::pow(start[i] - x[i] / A[i][i], 2);
				start[i] = x[i] / A[i][i];
				x[i] = b[i];
			}
		} while (std::sqrt(nrm) > eps && nrm < 1e10);

		std::cout << lp << "\n";
		return start;
	
	}

	template<class T, typename tmp_T = double>
	dynamic_vector<T>
		relax_iteration(
			const dynamic_matrix<T>&   A,
			const dynamic_vector<T>&   b,
			tmp_T                      eps, 
			dynamic_vector<T>          start = dynamic_vector<T>(),
			tmp_T                      omega = 0.5) {

		if (A.rows() != A.columns()) throw(std::invalid_argument("Sizes of matrix rows and columns must be equal"));
		if (A.rows() != b.size()) throw(std::invalid_argument("Matrix and vector have different sizes"));

		if (start.size() == 0) start.copy(b);
		else if (start.size() != b.size()) throw(std::invalid_argument("Vector and vector have different sizes"));

		size_t n = b.size();

		dynamic_vector<T> x(b);
		tmp_T nrm = 0;
		size_t lp = 0;
		do {
			++lp;
			nrm = 0;

			for (size_t j = 0; j != n; ++j){
				for (size_t i = j + 1; i < n; ++i) {
					x[j] -= A[j][i] * start[i];
				}

				for (size_t i = 0; i < j ; ++i) {
					x[j] -= A[j][i] * x[i];
				}

				x[j] *= omega;
				x[j] /= A[j][j];

				x[j] += (1 - omega) * start[j];
			}

			for (size_t i = 0; i != n; ++i) {
				nrm = std::pow(start[i] - x[i], 2);
				start[i] = x[i];
				x[i] = b[i];
			}
		} while (std::sqrt(nrm) > eps && nrm < 1e10);

		std::cout << lp << "\n";
		return start;

	}

	template<class T, typename tmp_T = double>
	dynamic_matrix<T> minus_E(dynamic_matrix<T>& A) {
		dynamic_matrix<T> A_new(A);
		size_t n = A.rows();
		for (size_t i = 0; i < n; ++i) {
			A_new[i][i] = A[i][i] - 1;
		}

		return A_new;
	}

	template<class T, typename tmp_T = double>
	dynamic_vector<T> simple_iter(
		const dynamic_matrix<T>&       A,
		const dynamic_vector<T>&       b,
		const tmp_T					   eps,
		const tmp_T					       tau = 0.05) {
		std::cout << "Test A: \n";
		mymath::utilities::print(A);
		//start here
		dynamic_matrix<T> E(A);
		dynamic_matrix<T>::diag(E, 1);

		std::cout << "Test E: \n";
		mymath::utilities::print(E);

		dynamic_vector<T> x(b);
		dynamic_vector<T> tmp_x(b);
		dynamic_vector<T> x_new(0, b.size());

		dynamic_matrix<T> C = E - tau * A;
		std::cout << "Test C: \n";
		mymath::utilities::print(C);

		//find tau cycle
		/*
		tmp_T tau1 = -1;
		tmp_T minim = 100;
		while (cube_norm(C) >= 1) {
			C = E - tau1 * A;
			minim = min(cube_norm(C), minim);
			tau1 += 0.00001;
			if (tau1 > 1) {std::cout << minim << "\n"; break;}
		}
		tau = tau1;
		*/

		std::cout << "norm C: " << cube_norm(C) << "\n";

		size_t k = 0;
		if (cube_norm(C) < 1) {
			while (cube_norm(x_new - x) > ((1 - cube_norm(C)) / cube_norm(C)) * eps) {
				tmp_x.move(x_new);
				x = multiply(C, x);
				x_new = tau * b + x;
				x.move(tmp_x);
				k += 1;
			}
			//std::cout << "vector x^ \n";
			//utilities::print(x_new);
		}
		else { std::cout << "||C|| >= 1\n\n"; }
		std::cout << k << "\n";
		//end here
		return x;
	}

#endif
}