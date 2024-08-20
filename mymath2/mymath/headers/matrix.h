#ifndef MYMATH_MATRIX
#define MYMATH_MATRIX
#include <cmath>

namespace mymath{


	// ------------------------------------------------------------------>	 CLASSES

	template<class T, size_t n, size_t m>
	class matrix{ 
	

		class row {

		};

		class const_row {

		};

	public:
		T values[n][m];

		static matrix<T, n, m>& fill(matrix<T, n, m>&, const T& elem = 0);

		static matrix<T, n, m>& eve(matrix<T, n, m>&);

		matrix<T, n, m>& operator+=(const matrix<T, n, m>&);

		matrix<T, n, m>& operator-=(const matrix<T, n, m>&);

		matrix<T, n, m>& operator*=(const T&);

		matrix<T, n, m>& operator/=(const T&);
	};


	

	

}





#endif
