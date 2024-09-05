#pragma once
#include "../../headers/matrix.h"

// ------------------------------------------------------------------>	 EXTERNAL OPERATORS

	// ------------------ MATRIX PLUS ------------------
namespace mymath {
	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator+(const matrix<T, n, m>& A, const matrix<T, n, m>& B) {
		matrix<T, n, m> tmp(A);
		tmp += B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator+(matrix<T, n, m>&& A, const matrix<T, n, m>& B) {
		A += B;
		return A;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator+(const matrix<T, n, m>& A, matrix<T, n, m>&& B) {
		B += A;
		return B;
	}

	// ------------------ MATRIX MINUS ------------------

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator-(const matrix<T, n, m>& A, const matrix<T, n, m>& B) {
		matrix<T, n, m> tmp(A);

		tmp -= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator-(matrix<T, n, m>&& A, const matrix<T, n, m>& B) {
		A -= B;
		return A;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator-(const matrix<T, n, m>& A, matrix<T, n, m>&& B) {
		B -= A;
		B *= -1;
		return B;
	}

	// ------------------ MATRIX NUM MUL ------------------

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator*(const matrix<T, n, m>& A, const T& B) {
		matrix<T, n, m> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator*(const T& B, const matrix<T, n, m>& A) {
		matrix<T, n, m> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator*(const matrix<T, n, m>&& A, const T& B) {
		A *= B;

		return A;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator*(const T& B, const matrix<T, n, m>&& A) {
		A *= B;

		return A;
	}


	// ------------------ MATRIX NUM DIV ------------------

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator/(const matrix<T, n, m>& A, const T& B) {
		matrix<T, n, m> tmp(A);

		tmp /= B;

		return tmp;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m> operator/(const matrix<T, n, m>&& A, const T& B) {
		A /= B;

		return A;
	}

}