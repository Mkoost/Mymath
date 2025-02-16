#pragma once

#include "../matrix.inl.h"
#include "../vector.inl.h"
#include "../dynamic_matrix.inl.h"
#include "../dynamic_vector.inl.h"
#include "../../details/__expr.inl"
#include "../../settings.h"
#include <cmath>


namespace mymath {



	namespace {
		template<class T, class U>
		void __make_jacobi(dynamic_matrix<T>& mat, dynamic_vector<T>& p, T t, const dynamic_vector<T>& x, U f, double eps) {
			size_t n = f.size();
			f[0];
			for (size_t i = 0; i < n; ++i) {
				p[i] = f[i](t, x);
			}

			dynamic_vector<T> tmp = x;

			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j) {
					tmp[j] = x[j] + eps;
					mat[i][j] = (f[i](t, tmp) - p[i]) / eps;
					tmp[j] = x[j];
				}

		}
	}

	template<class T, class U>
	dynamic_vector<T> eq_system_solve(T t, const dynamic_vector<T>& loc_p, U& f, double eps = 1e-8, size_t max_iter = 2) {

		dynamic_vector<T> x = loc_p;
		dynamic_vector<T> tmpp = loc_p;
		dynamic_matrix<T> jacobi(0, f.size(), f.size());
		dynamic_vector<T> v(0, f.size());
		double nrm = 1000;

		size_t iter = 0;

		do {
			++iter;
			__make_jacobi(jacobi, v, t, x, f, eps / 10);

			tmpp = x;
			v *= -1;
			v += multiply(jacobi, x);

			x = relax_iteration(jacobi, v, std::max(1e-11, eps / 10), 1.0, x);
			tmpp -= x;
			nrm = mymath::cube_norm(tmpp);
		} while (eps < nrm && iter <= max_iter);
		return x;
	}



}