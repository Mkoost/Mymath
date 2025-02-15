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
	dynamic_vector<T> zeidel_eq_system_solve(dynamic_vector<T> loc_p,  dynamic_matrix<T> localization, dynamic_vector<U> f, dynamic_vector<T> st_p, double eps, size_t max_iter=50) {
		
	

		for (size_t i = 0; i < max_iter; ++i) {
			
			return {0};
		}


	}



}