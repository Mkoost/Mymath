#pragma once
#include "../../headers/dynamic_matrix.h"

// ------------------------------------------------------------------>	 EXTERNAL OPERATORS

	// ------------------ MATRIX PLUS ------------------
namespace mymath {
	

	template<class T>
	dynamic_matrix<T> operator+(const dynamic_matrix<T>& A, const dynamic_matrix<T>& B) {
		dynamic_matrix<T> tmp(A);
		tmp += B;

		return tmp;
	};


	template<class T>
	dynamic_matrix<T> operator+(dynamic_matrix<T>&& A, const dynamic_matrix<T>& B) {
		A += B;
		return A;
	}


	template<class T>
	dynamic_matrix<T> operator+(const dynamic_matrix<T>& A, dynamic_matrix<T>&& B) {
		B += A;
		return B;
	}

	// ------------------ MATRIX MINUS ------------------

	template<class T>
	dynamic_matrix<T> operator-(const dynamic_matrix<T>& A, const dynamic_matrix<T>& B) {
		dynamic_matrix<T> tmp(A);

		tmp -= B;

		return tmp;
	}

	template<class T>
	dynamic_matrix<T> operator-(dynamic_matrix<T>&& A, const dynamic_matrix<T>& B) {
		A -= B;
		return A;
	}

	template<class T>
	dynamic_matrix<T> operator-(const dynamic_matrix<T>& A, dynamic_matrix<T>&& B) {
		B -= A;
		B *= -1;
		return B;
	}

	// ------------------ MATRIX NUM MUL ------------------

	template<class T>
	constexpr dynamic_matrix<T> operator*(const dynamic_matrix<T>& A, const T& B) {
		dynamic_matrix<T> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T>
	constexpr dynamic_matrix<T> operator*(const T& B, const dynamic_matrix<T>& A) {
		dynamic_matrix<T> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T>
	constexpr dynamic_matrix<T> operator*(const dynamic_matrix<T>&& A, const T& B) {
		A *= B;

		return A;
	}

	template<class T>
	constexpr dynamic_matrix<T> operator*(const T& B, const dynamic_matrix<T>&& A) {
		A *= B;

		return A;
	}


	// ------------------ MATRIX NUM DIV ------------------

	template<class T>
	constexpr dynamic_matrix<T> operator/(const dynamic_matrix<T>& A, const T& B) {
		dynamic_matrix<T> tmp(A);

		tmp /= B;

		return tmp;
	}

	template<class T>
	constexpr dynamic_matrix<T> operator/(const dynamic_matrix<T>&& A, const T& B) {
		A /= B;

		return A;
	}


}