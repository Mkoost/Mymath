#pragma once
#include "matrix.inl.h"
#include "vector.inl.h"
namespace mymath {
	namespace utilities {
#if defined(MYMATH_QUATERNION_STATE)
		template<class T, size_t n, size_t m>
		void print(const mymath::matrix<T, n, m>& A) {
			for (int i = 0; i != n; ++i) {
				for (int j = 0; j != m; ++j)
					std::cout << A[i][j] << " ";
				std::cout << "\n";
			}
			std::cout << "\n";
		}
#endif

#if defined(MYMATH_MATRIX_STATE)
		template<class T>
		void print(const mymath::quaternion<T>& q) {
			std::cout << q.w << " " << q.x << " " << q.y << " " << q.z;
		}
#endif

#if defined(MYMATH_VECTOR_STATE)
		template<class T, size_t n>
		void print(const mymath::vector<T, n>& vec) {
			for (int i = 0; i != n; ++i)
				std::cout << vec[i] << " ";
			std::cout << "\n";
		}
#endif

	}
}
