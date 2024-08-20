
#pragma once


//#include "../headers/matrix.h"

namespace mymath {


	// ------------------------------------------------------------------>	 CLASS METHODS 

	// ------------------ STATIC METHODS ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m>& matrix<T, n, m>::fill(matrix<T, n, m>& mat, const T& some) {
		T* ptr = reinterpret_cast<T*>(&mat.values);
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			ptr[i] = some;
		return mat;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m>& matrix<T, n, m>::eve(matrix<T, n, m>& mat) {
		fill(mat, 0);
		for (size_t i = 0; i != std::min(n, m); ++i) {
			mat.values[i][i] = 1;
		}
		return mat;
	}

	// ------------------ OPERATORS ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m>& matrix<T, n, m>::operator+=(const matrix<T, n, m>& mat) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] += reinterpret_cast<const T*>(&mat.values)[i];
		return *this;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m>& matrix<T, n, m>::operator-=(const matrix<T, n, m>& mat) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] -= reinterpret_cast<const T*>(&mat.values)[i];
		return *this;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m>& matrix<T, n, m>::operator*=(const T& a) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] *= a;
		return *this;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m>& matrix<T, n, m>::operator/=(const T& a) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] /= a;
		return *this;
	}



	// ------------------------------------------------------------------>	 EXTERNAL OPERATORS

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
	matrix<T, n, m> operator+(const matrix<T, n, m>& A, matrix<T, n, m>&& B) {
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