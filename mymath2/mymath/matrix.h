//#if defined(MYMATH_QUATERNION_STATE) && (MYMATH_MATRIX_STATE == 1)
//#endif

#ifndef MYMATH_MATRIX
#define MYMATH_MATRIX
#include <iostream>

namespace mymath{


	// ------------------------------------ CLASSES ------------------------------------

	template<class T, size_t n, size_t m>
	class matrix{ 
	

		class row {

		};

		class const_row {

		};

	public:
		T values[n][m];

		static matrix<T,  n,  m>& fill(matrix<T, n, m>& mat, const T& some = 0) {
			T* ptr = reinterpret_cast<T*>(&mat.values);
			size_t k = n * m;
			for (size_t i = 0; i != k; ++i)
				ptr[i] = some;
			return mat;
		}

		static matrix<T,  n,  m>& eve(matrix<T, n, m>& mat) {
			zeros(mat);
			for (size_t i = 0; i != n; ++i){
				mat.values[i][i] = 1;
			}
			return mat;
		}

		matrix<T, n, m>& operator+=(const matrix<T, n, m>& mat) {
			size_t k = n * m;
			for (size_t i = 0; i != k; ++i)
				reinterpret_cast<T*>(values)[i] += reinterpret_cast<const T*>(&mat.values)[i];
			return *this;
		}

		matrix<T, n, m>& operator-=(const matrix<T, n, m>& mat) {
			size_t k = n * m;
			for (size_t i = 0; i != k; ++i)
				reinterpret_cast<T*>(values)[i] -= reinterpret_cast<const T*>(&mat.values)[i];
			return *this;
		}

		matrix<T, n, m>& operator*=(const T& a) {
			size_t k = n * m;
			for (size_t i = 0; i != k; ++i)
				reinterpret_cast<T*>(values)[i] *= a;
			return *this;
		}

		matrix<T, n, m>& operator/=(const T& a) {
			size_t k = n * m;
			for (size_t i = 0; i != k; ++i)
				reinterpret_cast<T*>(values)[i] /= a;
			return *this;
		}


	};


	// ------------------------------------ OPERATORS ------------------------------------

	// ------------------ MATRIX PLUS ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator+(const matrix<T, n, m>& A, const matrix<T, n, m>& B) {
		matrix<T, n, m> tmp(A);
		tmp += B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator+(matrix<T, n, m>&& A, const matrix<T, n, m>& B) {
		A += B;
		return A;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator+(const matrix<T, n, m>& A,  matrix<T, n, m>&& B) {
		B += A;
		return B;
	}

	// ------------------ MATRIX MINUS ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator-(const matrix<T, n, m>& A, const matrix<T, n, m>& B) {
		matrix<T, n, m> tmp(A);

		tmp -= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator-(matrix<T, n, m>&& A, const matrix<T, n, m>& B) {
		A -= B;
		return A;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator-(const matrix<T, n, m>& A, matrix<T, n, m>&& B) {
		B -= A;
		B *= -1;
		return B;
	}

	// ------------------ MATRIX NUM MUL ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator*(const matrix<T, n, m>& A, const T& B) {
		matrix<T, n, m> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator*(const T& B, const matrix<T, n, m>& A) {
		matrix<T, n, m> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator*(const matrix<T, n, m>&& A, const T& B) {
		A *= B;

		return A;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator*(const T& B, const matrix<T, n, m>&& A) {
		A *= B;

		return A;
	}

	// ------------------ MATRIX NUM DIV ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator/(const matrix<T, n, m>& A, const T& B) {
		matrix<T, n, m> tmp(A);

		tmp /= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m> operator/(const matrix<T, n, m>&& A, const T& B) {
		A /= B;

		return A;
	}

	

}





#endif
