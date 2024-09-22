#pragma once
#include "./../../headers/dynamic_matrix.h"
#include "./../../details/__expr.inl"

namespace mymath {
	// ------------------------------------------------------------------>	 CLASS METHODS 

	// ------------------ STATIC METHODS ------------------

	template<class T>
	const dynamic_matrix<T>& dynamic_matrix<T>::fill(const dynamic_matrix<T>& mat, const T& some) {
		T* ptr = reinterpret_cast<T*>(const_cast<dynamic_matrix<T>*>(&mat)->values);
		for (size_t i = 0, k = mat.size(); i != k; ++i)
			ptr[i] = some;
		return mat;
	}

	template<class T>
	const dynamic_matrix<T>& dynamic_matrix<T>::diag(const dynamic_matrix<T>& mat, const T& some) {
		fill(mat, 0);
		T* ptr = reinterpret_cast<T*>(const_cast<dynamic_matrix<T>*>(&mat)->values);
		for (size_t i = 0; i != min(mat.columns(), mat.rows()); ++i) {
			ptr[i * (mat.columns() + 1)] = some;
		}
		return mat;
	}
	
	// ------------------ CONSTRUCTOR / DESTRUCTOR ------------------
	

	template<class T>
	dynamic_matrix<T>::dynamic_matrix(const T& some, size_t n, size_t m){
		n_ = n;
		m_ = m;
		values = new T[n_ * m_];
		fill(*this, some);
	}

	template<class T>
	dynamic_matrix<T>::dynamic_matrix(const std::initializer_list<std::initializer_list<T>>& list) {
		
		n_ = list.size();
		for (auto i : list)
			m_ = max(i.size(), m_);

		values = new T[n_ * m_];

		size_t n = 0, m = 0;
		for (auto i : list) {
			for (auto j : i)
				values[(m++) + m_ * n] = j;
			while (m <= m_ - 1)
				values[(m++) + m_ * n] = 0;
			m = 0;
			++n;
		}
	}

	template<class T>
	template<class U>
	dynamic_matrix<T>::dynamic_matrix(const dynamic_matrix<U>& A) {
		copy(A);
	}

	template<class T>
	dynamic_matrix<T>::dynamic_matrix(const dynamic_matrix<T>& A) {
		copy(A);
	}

	template<class T>
	dynamic_matrix<T>::dynamic_matrix(T* A, size_t n, size_t m) {
		values = A;
		n_ = n; m_ = m;
	}

	template<class T>
	dynamic_matrix<T>::dynamic_matrix(dynamic_matrix<T>&& A) {
		move(A);
	}

	template<class T>
	dynamic_matrix<T>::~dynamic_matrix() {
		if (values) delete[] values;
	}

	// ------------------ METHODS ------------------

	template<class T>
	size_t dynamic_matrix<T>::columns() const { return m_; }

	template<class T>
	size_t dynamic_matrix<T>::rows() const { return n_; }

	template<class T>
	size_t dynamic_matrix<T>::size() const { return n_ * m_; }

	template<class T>
	unsigned dynamic_matrix<T>::id() const { return 1u; }


	template<class T>
	template<class U>
	dynamic_matrix<T>& dynamic_matrix<T>::copy(const dynamic_matrix<U>& A) {
		if (values) delete[] values;

		n_ = A.rows();
		m_ = A.columns();

		values = new T[A.size()];

		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] = static_cast<T>(A.values[i]);

		return *this;
	};

	template<class T>
	dynamic_matrix<T>& dynamic_matrix<T>::move(dynamic_matrix<T>& A) noexcept {
		if (values) delete[] values;

		values = A.values;
		n_ = A.rows();
		m_ = A.columns();

		A.values = nullptr;
		A.m_ = A.n_ = 0;
		return *this;
	}


	/*template<class T>
	constexpr typename dynamic_matrix<T>::iterator dynamic_matrix<T>::begin() noexcept {
		return iterator(&values);
	};*/

	/*template<class T>
	constexpr typename dynamic_matrix<T>::iterator dynamic_matrix<T>::end() noexcept {
		return iterator(& values[n - 1][m - 1] + 1);
	};*/

	/*template<class T>
	constexpr typename dynamic_matrix<T>::const_iterator dynamic_matrix<T>::begin() const noexcept  {
		return const_iterator(&values);
	};*/

	/*template<class T>
	constexpr typename dynamic_matrix<T>::const_iterator dynamic_matrix<T>::end() const noexcept {
		return const_iterator(& values[n - 1][m - 1] + 1);
	};*/

	// ------------------ OPERATORS ------------------

	template<class T>
	constexpr dynamic_matrix<T>& dynamic_matrix<T>::operator+=(const dynamic_matrix<T>& mat) { 
		if (mat.m_ != m_ && mat.n_ != m_) throw(std::invalid_argument("different sizes of matrices"));
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] += mat.values[i];
		return *this;
	}

	template<class T>
	constexpr dynamic_matrix<T>& dynamic_matrix<T>::operator-=(const dynamic_matrix<T>& mat) {
		if (mat.m_ != m_ && mat.n_ != m_) throw(std::invalid_argument("different sizes of matrices"));
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] -= mat.values[i];
		return *this;
	}

	template<class T>
	constexpr dynamic_matrix<T>& dynamic_matrix<T>::operator*=(const T& a) {
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] *= a;
		return *this;
	}

	template<class T>
	constexpr dynamic_matrix<T>& dynamic_matrix<T>::operator/=(const T& a) {
		for (size_t i = 0, k = size(); i != k; ++i)
			values[i] /= a;
		return *this;
	}

	template<class T>
	constexpr dynamic_matrix<T>& dynamic_matrix<T>::operator=(const dynamic_matrix<T>&A) {
		copy(A);
		return *this;
	}

	template<class T>
	constexpr dynamic_matrix<T>& dynamic_matrix<T>::operator=( dynamic_matrix<T>&& A) {
		move(A);
		return *this;
	}

	template<class T>
	constexpr typename dynamic_matrix<T>::row dynamic_matrix<T>::operator[](size_t i) {
		row r(values + i * m_);
		return r;
	};

	template<class T>
	constexpr typename dynamic_matrix<T>::const_row dynamic_matrix<T>::operator[](size_t i) const {
		const_row r(values + i * m_);
		return r;
	};

}