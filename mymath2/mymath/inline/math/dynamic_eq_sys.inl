#pragma once

#include "../matrix.inl.h"
#include "../vector.inl.h"
#include "../dynamic_matrix.inl.h"
#include "../dynamic_vector.inl.h"
#include "../../details/__expr.inl"
#include "../../settings.h"
#include <cmath>


namespace mymath {

	//                ^ v_3
	//                |
	//                |
	//                |
	//                |
	//                |
	//             P  0----------------> v_2
	//               /
	//              /
	//             /
	//            v
	//            v_1
	//                                                       
	// localization matrix: (v_1, v_2, v_3), where is v = (x1, x2, x3)^T
	// localization point: is P 
	//
	// D1 Example:
	// 
	//          y
	//			^
	//          |                              /
	//          |                             |          
	//          |                            /
	//          |                          /
	//          |                        --
	//          |				        /
	//          |                     --
	//          |                    /
	//          |                 ---
	//          |             ----
	//          |   P        /           v_1
	//          0---0-------0------------>----------> x
	//          |  x_1     /            x_2
	//          |      ----
	//          |   ---   
	//          |  /
	//          |--
	//        / |
	//      --   |
	//     /     |
	//
	//	P = x_1
	//	v1 = x_2 - x_1
	//
	// D2 Example:
	//
	//          y
	//			^
	//          |                    v_2        /
	//          |                     ^        |          
	//          |                    /        /
	//          |                   /       /
	//          |         ---------/-----0---------
	//          |				  /     /
	//     y_1  0             p  0    --
	//          |                 \  /
	//          |                 -\-
	//          |             ----  \
	//          |                    > v_1  
	//          0---------------0------------------> x
	//                         x_1
	// P = (x_1, y_1)
	// Loc = (v_1, v_2)       


		template<class T, class U>
		void make_jacobi(dynamic_matrix<T>& mat, dynamic_vector<T>& p, T t, const dynamic_vector<T>& x, const dynamic_vector<U>& fl, const dynamic_vector<U>& fr, double eps) {
			size_t n = fl.size();
			for (size_t i = 0; i < n; ++i){
				p[i] = fl[i](t, x) - fr[i](t, x);
			}

			dynamic_vector<T> tmp = x;

			for(size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j) {
					tmp[j] = x[j] + eps;
					mat[i][j] = (fl[i](t, tmp) - fr[i](t, tmp) - p[i]) / eps;
					tmp[j] = x[j];
				}

		}


	template<class T, class U>
	dynamic_vector<T> eq_system_solve(T t, const dynamic_vector<T>& loc_p, const dynamic_vector<U>& fl, const dynamic_vector<U>& fr, double eps = 1e-8, size_t max_iter = 2) {
		
		if ((fl.size() != fr.size()))
			throw(std::invalid_argument("Number of left and right equations are not equal"));

		dynamic_vector<T> x = loc_p;
		dynamic_matrix<T> jacobi(0, fl.size(), fl.size());
		dynamic_vector<T> v(0, fl.size());
		double nrm = 1000;
		size_t iter = 0;

		do {
			++iter;
			make_jacobi(jacobi, v, t, x, fl, fr, eps / 10);
			
			
			v *= -1;
			v += multiply(jacobi, x);

			x = relax_iteration(jacobi, v, std::max( 1e-11, eps / 10), 1.0, x);
			

		}while(eps < nrm && iter <= max_iter);
		return x;
	}



}