#pragma once
#include "./../../headers/dynamic_vector.h"
#include "./../../details/__expr.inl"
#include <stdexcept>


namespace mymath {
	// ------------------------------------------------------------------>	 CLASS METHODS 

	// ------------------ STATIC METHODS ------------------

	template<class T>
	const dynamic_vector<T>& dynamic_vector<T>::fill(const dynamic_vector<T>& mat, const T& some) {
		T* ptr = reinterpret_cast<T*>(const_cast<dynamic_vector<T>*>(&mat)->values);
		for (size_t i = 0, k = mat.size(); i != k; ++i)
			ptr[i] = some;
		return mat;
	}

	
	// ------------------ CONSTRUCTOR / DESTRUCTOR ------------------
	

	template<class T>
	dynamic_vector<T>::dynamic_vector(const T& some, size_t n){
		n_ = n;
		values = new T[n_];
		fill(*this, some);
	}

	template<class T>
	dynamic_vector<T>::dynamic_vector(const std::initializer_list<T>& list) {
		
		n_ = list.size();

		values = new T[size()];

		size_t j = 0;
		for (auto i : list) 
			values[j++] = i;
	}

	template<class T>
	template<class U>
	dynamic_vector<T>::dynamic_vector(const dynamic_vector<U>& A) {
		copy(A);
	}

	template<class T>
	dynamic_vector<T>::dynamic_vector(const dynamic_vector<T>& A) {
		copy(A);
	}

	template<class T>
	dynamic_vector<T>& dynamic_vector<T>::move(T* A, size_t n) noexcept {
		values = A;
		n_ = n;
		return *this;
	}

	template<class T>
	dynamic_vector<T>::dynamic_vector(dynamic_vector<T>&& A) {
		move(A);
	}

	template<class T>
	dynamic_vector<T>::~dynamic_vector() {
		if (values) delete[] values;
	}

	// ------------------ METHODS ------------------

	template<class T>
	size_t dynamic_vector<T>::size() const { return n_; }

	template<class T>
	unsigned dynamic_vector<T>::id() const { return 2u; }


	template<class T>
	template<class U>
	dynamic_vector<T>& dynamic_vector<T>::copy(const dynamic_vector<U>& A) {
		if (values) delete[] values;

		n_ = A.size();

		values = new T[A.size()];

		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] = static_cast<T>(A.values[i]);

		return *this;
	};

	template<class T>
	dynamic_vector<T>& dynamic_vector<T>::move(dynamic_vector<T>& A) noexcept {
		if (values) delete[] values;

		values = A.values;
		n_ = A.size();

		A.values = nullptr;
		A.n_ = 0;
		return *this;
	}


	/*template<class T>
	constexpr typename dynamic_vector<T>::iterator dynamic_vector<T>::begin() noexcept {
		return iterator(&values);
	};*/

	/*template<class T>
	constexpr typename dynamic_vector<T>::iterator dynamic_vector<T>::end() noexcept {
		return iterator(& values[n - 1][m - 1] + 1);
	};*/

	/*template<class T>
	constexpr typename dynamic_vector<T>::const_iterator dynamic_vector<T>::begin() const noexcept  {
		return const_iterator(&values);
	};*/

	/*template<class T>
	constexpr typename dynamic_vector<T>::const_iterator dynamic_vector<T>::end() const noexcept {
		return const_iterator(& values[n - 1][m - 1] + 1);
	};*/

	// ------------------ OPERATORS ------------------

	template<class T>
	constexpr dynamic_vector<T>& dynamic_vector<T>::operator+=(const dynamic_vector<T>& mat) { 
		if (mat.size() > size()) throw(std::invalid_argument("The argument's vector size is bigger than this size"));
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] += mat.values[i];
		return *this;
	}

	template<class T>
	constexpr dynamic_vector<T>& dynamic_vector<T>::operator-=(const dynamic_vector<T>& mat) {
		if (mat.size() > size()) throw(std::invalid_argument("The argument's vector size is bigger than this size"));
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] -= mat.values[i];
		return *this;
	}

	template<class T>
	constexpr dynamic_vector<T>& dynamic_vector<T>::operator*=(const T& a) {
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] *= a;
		return *this;
	}

	template<class T>
	constexpr dynamic_vector<T>& dynamic_vector<T>::operator/=(const T& a) {
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] /= a;
		return *this;
	}

	template<class T>
	constexpr dynamic_vector<T>& dynamic_vector<T>::operator=(const dynamic_vector<T>&A) {
		copy(A);
		return *this;
	}

	template<class T>
	constexpr dynamic_vector<T>& dynamic_vector<T>::operator=( dynamic_vector<T>&& A) {
		move(A);
		return *this;
	}

	template<class T>
	constexpr T& dynamic_vector<T>::operator[](size_t i) {
		return values[i];
	};

	template<class T>
	constexpr T dynamic_vector<T>::operator[](size_t i) const {
		return values[i];
	};

}