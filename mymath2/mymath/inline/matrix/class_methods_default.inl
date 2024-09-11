#pragma once
#include "./../../headers/matrix.h"
#include "./../../details/__expr.inl"

namespace mymath {
	// ------------------------------------------------------------------>	 CLASS METHODS 

	// ------------------ STATIC METHODS ------------------

	template<class T, size_t n, size_t m>
	const matrix<T, n, m>& matrix<T, n, m>::fill(const matrix<T, n, m>& mat, const T& some) {
		T* ptr = reinterpret_cast<T*>(const_cast<matrix<T, n, m> *>(&mat)->values);
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			ptr[i] = some;
		return mat;
	}

	template<class T, size_t n, size_t m>
	const matrix<T, n, m>& matrix<T, n, m>::diag(const matrix<T, n, m>& mat, const T& some) {
		fill(mat, 0);
		T* ptr = reinterpret_cast<T*>(const_cast<matrix<T, n, m> *>(&mat)->values);
		for (size_t i = 0; i != min(n, m); ++i) {
			ptr[i * (m + 1)] = some;
		}
		return mat;
	}
	// ------------------ METHODS ------------------

	template<class T, size_t n, size_t m>
	template<class U>
	matrix<T, n, m>& matrix<T, n, m>::copy(const matrix<U, n, m>& A) {
		T* ptr1 = (T*)values;
		U* ptr2 = (U*)A.values[0][0];
		for (size_t i = 0; i != n * m; ++i)
			ptr1[i] = ptr2[i];
		return *this;
	};

	template<class T, size_t n, size_t m>
	constexpr typename matrix<T, n, m>::iterator matrix<T, n, m>::begin() noexcept {
		return iterator(&values);
	};

	template<class T, size_t n, size_t m>
	constexpr typename matrix<T, n, m>::iterator matrix<T, n, m>::end() noexcept {
		return iterator(& values[n - 1][m - 1] + 1);
	};

	template<class T, size_t n, size_t m>
	constexpr typename matrix<T, n, m>::const_iterator matrix<T, n, m>::begin() const noexcept  {
		return const_iterator(&values);
	};

	template<class T, size_t n, size_t m>
	constexpr typename matrix<T, n, m>::const_iterator matrix<T, n, m>::end() const noexcept {
		return const_iterator(& values[n - 1][m - 1] + 1);
	};
	// ------------------ CONSTRUCTOR ------------------

	template<class T, size_t n, size_t m>
	matrix<T, n, m>::matrix(const T& some){
		fill(*this, some);
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m>::matrix(const std::initializer_list<T>& list) {
		size_t k = 0;
		T* ptr = &values[0][0];
		for (auto i : list) 
			if (k != n * m){
				ptr[k++] = i;
			}

		while (k != n * m)
			ptr[k++] = 0;
	}

	template<class T, size_t n, size_t m>
	matrix<T, n, m>::matrix(const std::initializer_list<std::initializer_list<T>>& list) {
		size_t k = 0;
		T* ptr = &values[0][0];
		for (auto i : list)
			for (auto j : i)
			if (k != n * m) {
				ptr[k++] = j;
			}

		while (k != n * m)
			ptr[k++] = 0;
	}

	template<class T, size_t n, size_t m>
	template<class U>
	matrix<T, n, m>::matrix(const matrix<U, n, m>& A) {
		copy(A);
	}
	// ------------------ OPERATORS ------------------

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m>& matrix<T, n, m>::operator+=(const matrix<T, n, m>& mat) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] += reinterpret_cast<const T*>(&mat.values)[i];
		return *this;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m>& matrix<T, n, m>::operator-=(const matrix<T, n, m>& mat) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] -= reinterpret_cast<const T*>(&mat.values)[i];
		return *this;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m>& matrix<T, n, m>::operator*=(const T& a) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] *= a;
		return *this;
	}

	template<class T, size_t n, size_t m>
	constexpr matrix<T, n, m>& matrix<T, n, m>::operator/=(const T& a) {
		size_t k = n * m;
		for (size_t i = 0; i != k; ++i)
			reinterpret_cast<T*>(values)[i] /= a;
		return *this;
	}

	template<class T, size_t n, size_t m>
	constexpr typename matrix<T, n, m>::row matrix<T, n, m>::operator[](size_t i) {
		row r;
		r.ptr = &values[i][0];
		return r;
	};

	template<class T, size_t n, size_t m>
	constexpr typename matrix<T, n, m>::const_row matrix<T, n, m>::operator[](size_t i) const {
		const_row r;
		r.ptr = &values[i][0];
		return r;
	};

}