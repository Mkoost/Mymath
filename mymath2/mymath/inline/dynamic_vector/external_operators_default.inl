#pragma once
#include "../../headers/dynamic_vector.h"

// ------------------------------------------------------------------>	 EXTERNAL OPERATORS

	// ------------------ VECTOR PLUS ------------------
namespace mymath {
	

	template<class T>
	dynamic_vector<T> operator+(const dynamic_vector<T>& A, const dynamic_vector<T>& B) {
		dynamic_vector<T> tmp(A);
		tmp += B;

		return tmp;
	};


	template<class T>
	dynamic_vector<T> operator+(dynamic_vector<T>&& A, const dynamic_vector<T>& B) {
		A += B;
		return A;
	}


	template<class T>
	dynamic_vector<T> operator+(const dynamic_vector<T>& A, dynamic_vector<T>&& B) {
		B += A;
		return B;
	}

	// ------------------ VECTOR MINUS ------------------

	template<class T>
	dynamic_vector<T> operator-(const dynamic_vector<T>& A, const dynamic_vector<T>& B) {
		dynamic_vector<T> tmp(A);

		tmp -= B;

		return tmp;
	}

	template<class T>
	dynamic_vector<T> operator-(dynamic_vector<T>&& A, const dynamic_vector<T>& B) {
		A -= B;
		return A;
	}

	template<class T>
	dynamic_vector<T> operator-(const dynamic_vector<T>& A, dynamic_vector<T>&& B) {
		B -= A;
		B *= -1;
		return B;
	}

	// ------------------ VECTOR NUM MUL ------------------

	template<class T>
	constexpr dynamic_vector<T> operator*(const dynamic_vector<T>& A, const T& B) {
		dynamic_vector<T> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T>
	constexpr dynamic_vector<T> operator*(const T& B, const dynamic_vector<T>& A) {
		dynamic_vector<T> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T>
	constexpr dynamic_vector<T> operator*(const dynamic_vector<T>&& A, const T& B) {
		A *= B;

		return A;
	}

	template<class T>
	dynamic_vector<T> operator*(const T& B, dynamic_vector<T>&& A) {
		A *= B;

		return A;
	}


	// ------------------ VECTOR NUM DIV ------------------

	template<class T>
	constexpr dynamic_vector<T> operator/(const dynamic_vector<T>& A, const T& B) {
		dynamic_vector<T> tmp(A);

		tmp /= B;

		return tmp;
	}

	template<class T>
	constexpr dynamic_vector<T> operator/(const dynamic_vector<T>&& A, const T& B) {
		A /= B;

		return A;
	}

	// ------------------ VECTOR MUL VECTOR ------------------

	template<class T>
	constexpr auto operator*(const dynamic_vector<T>& A, const dynamic_vector<T>& B) -> decltype(A[0] * B[0]) {
		auto tmp = A[0] * B[0];
		for (size_t i = 0, k = min(A.size(), B.size()); i != k; ++i)
			tmp += A[i] * B[i];
		return tmp;
	}

}