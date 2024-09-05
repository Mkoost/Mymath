#pragma once
#include "../../headers/vector.h"

namespace mymath {
// ------------------------------------------------------------------>	 EXTERNAL OPERATORS

	// ------------------ VECTOR PLUS ------------------

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator+(const  vector<T, n>& A, const  vector<U, n>& B) {
		vector<T, n> tmp(A);
		tmp += B;

		return tmp;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator+(vector<T, n>&& A, const  vector<U, n>& B) {
		A += B;
		return A;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator+(const  vector<T, n>& A, vector<U, n>&& B) {
		B += A;
		return B;
	}

	// ------------------ VECTOR MINUS ------------------
	template<class T, class U, size_t n>
	constexpr vector<T, n> operator-(const  vector<T, n>& A, const  vector<U, n>& B) {
		vector<T, n> tmp(A);
		tmp -= B;

		return tmp;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator-(vector<T, n>&& A, const  vector<U, n>& B) {
		A -= B;
		return A;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator-(const  vector<U, n>& A, vector<T, n>&& B) {
		B -= A;
		return B;
	}

	// ------------------ VECTOR NUM MUL ------------------

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator*(const vector<T, n>& A, const U& B) {
		vector<T, n> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator*(const U& B, const vector<T, n>& A) {
		vector<T, n> tmp(A);

		tmp *= B;

		return tmp;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator*(const vector<T, n>&& A, const U& B) {
		A *= B;

		return A;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator*(const U& B, const vector<T, n>&& A) {
		A *= B;

		return A;
	}
	

	// ------------------ VECTOR NUM DIV ------------------

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator/(const vector<T, n>& A, const U& B) {
		vector<T, n> tmp(A);

		tmp /= B;

		return tmp;
	}

	template<class T, class U, size_t n>
	constexpr vector<T, n> operator/(const vector<T, n>&& A, const U& B) {
		A /= B;

		return A;
	}

	// ------------------ VECTOR MUL VECTOR ------------------
	template<class T, class U, size_t n>
	constexpr auto operator*(const vector<T, n>& A, const vector<U, n>& B) -> decltype(A[0] * B[0]) {
		auto tmp = A[0] * B[0];
		for (size_t i = 1; i != n; ++i)
			tmp += A[i] * B[i];
		return tmp;
	}
}