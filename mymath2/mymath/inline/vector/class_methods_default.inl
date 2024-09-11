#pragma once
#include "./../../headers/vector.h"
#include "./../../details/__expr.inl"

namespace mymath {
	// ------------------------------------------------------------------>	 CLASS METHODS 

	// ------------------ STATIC METHODS ------------------

	template<class T, size_t n>
	const vector<T, n>& vector<T, n>::fill(const vector<T, n>& vec, const T& some) {
		T* ptr = reinterpret_cast<T*>(const_cast<vector<T, n> *>(&vec)->values);
		size_t k = n;
		for (size_t i = 0; i != k; ++i)
			ptr[i] = some;
		return vec;
	}

	// ------------------ METHODS ------------------
	template<class T, size_t n>
	template<class U>
	vector<T, n>& vector<T, n>::copy(const vector<U, n>& A) {
		for (size_t i = 0; i != n; ++i)
			values[i] = A.values[i];
		return *this;
	};

	template<class T, size_t n>
	constexpr typename vector<T, n>::iterator vector<T, n>::begin() noexcept {
		return iterator(&values);
	};

	template<class T, size_t n>
	constexpr typename vector<T, n>::iterator vector<T, n>::end() noexcept {
		return iterator(&values[n - 1] + 1);
	};

	template<class T, size_t n>
	constexpr typename vector<T, n>::const_iterator vector<T, n>::begin() const noexcept {
		return const_iterator(&values);
	};

	template<class T, size_t n>
	constexpr typename vector<T, n>::const_iterator vector<T, n>::end() const noexcept {
		return const_iterator(&values[n - 1]);
	};

	// ------------------ CONSTRUCTOR ------------------

	template<class T, size_t n>
	vector<T, n>::vector(const T& some) {
		fill(*this, some);
	}

	template<class T, size_t n>
	vector<T, n>::vector(const std::initializer_list<T>& list) {
		size_t k = 0;
		T* ptr = &values[0];
		for (auto i : list)
			if (k != n) {
				ptr[k++] = i;
			}

		while (k != n)
			ptr[k++] = 0;
	}

	template<class T, size_t n>
	template<class U>
	vector<T, n>::vector(const vector<U, n>& vec) {
		copy(vec);
	};

	// ------------------ OPERATORS ------------------

	template<class T, size_t n>
	template<class U>
	constexpr vector<T, n>& vector<T, n>::operator+=(const vector<U, n>& vec) {
		size_t k = n;
		for (size_t i = 0; i != k; ++i)
			values[i] += vec[i];
		return *this;
	}

	template<class T, size_t n>
	template<class U>
	constexpr vector<T, n>& vector<T, n>::operator-=(const vector<U, n>& vec) {
		size_t k = n;
		for (size_t i = 0; i != k; ++i)
			values[i] -= vec[i];
		return *this;
	}

	template<class T, size_t n>
	template<class U>
	constexpr vector<T, n>& vector<T, n>::operator*=(const U& a) {
		size_t k = n;
		for (size_t i = 0; i != k; ++i)
			values[i] *= a;
		return *this;
	}

	template<class T, size_t n>
	template<class U>
	constexpr vector<T, n>& vector<T, n>::operator/=(const U& a) {
		size_t k = n;
		for (size_t i = 0; i != k; ++i)
			values[i] /= a;
		return *this;
	}

	template<class T, size_t n>
	constexpr T& vector<T, n>::operator[](size_t i) {
		return values[i];
	};

	template<class T, size_t n>
	constexpr T vector<T, n>::operator[](size_t i) const {
		return values[i];
	};

}