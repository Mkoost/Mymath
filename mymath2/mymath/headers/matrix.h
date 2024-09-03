#pragma once


#define __GEN_MYMATH_MAT_PSEUDONYMS_10931(postfix, n, m) using dmat##postfix = matrix<double, n, m>; \
 using fmat##postfix = matrix<float, n, m>; \
using imat##postfix = matrix<int, n, m>; \
using lmat##postfix = matrix<long, n, m>; \
using llmat##postfix = matrix<long long, n, m>; \
using umat##postfix = matrix<unsigned, n, m>; \
using uimat##postfix = matrix<unsigned int, n, m>; \
using ulmat##postfix = matrix<unsigned long, n, m>; \
using ullmat##postfix = matrix<unsigned long long, n, m>;




namespace mymath{
	// ------------------------------------------------------------------>	 CLASSES

	template<class T, size_t n, size_t m>
	class matrix{ 
		
		
		struct row {
			T* ptr;
			T& operator[](size_t);
		};

		struct const_row {
			T const * ptr;
			T operator[](size_t) const;
		};

		T values[n][m];

	public:

		using iterator = typename T*;
		using const_iterator = typename T const *;

		matrix(const T& x);

		constexpr iterator begin() noexcept;
		constexpr iterator end() noexcept;
#if 0
		reverse_iterator rbegin();
		reverse_iterator rend();
#endif

		constexpr const_iterator  begin() const noexcept;
		constexpr const_iterator end() const noexcept;

#if 0
		reverse_const_iterator rbegin() const;
		reverse_const_iterator rend() const;
#endif

		/* WARNING: Can change const object */
		static const matrix<T, n, m>& fill(const matrix<T, n, m>&, const T& elem = 0);

		/* WARNING: Can change const object */
		static const matrix<T, n, m>& diag(const matrix<T, n, m>&, const T& elem = 1);

		constexpr matrix<T, n, m>& operator+=(const matrix<T, n, m>&);

		constexpr matrix<T, n, m>& operator-=(const matrix<T, n, m>&);

		constexpr matrix<T, n, m>& operator*=(const T&);

		constexpr matrix<T, n, m>& operator/=(const T&);

		constexpr row operator[](size_t);

		constexpr const_row operator[](size_t) const;
	};

	// ------------------------------------------------------------------>	 PSEUDONYMS

	__GEN_MYMATH_MAT_PSEUDONYMS_10931(2, 2, 2)
	__GEN_MYMATH_MAT_PSEUDONYMS_10931(3, 3, 3)
	__GEN_MYMATH_MAT_PSEUDONYMS_10931(4, 4, 4)
	__GEN_MYMATH_MAT_PSEUDONYMS_10931(5, 5, 5)

}

#undef __GEN_MYMATH_MAT_PSEUDONYMS_10931