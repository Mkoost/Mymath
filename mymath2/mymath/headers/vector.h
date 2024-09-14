#pragma once
#include <initializer_list>

namespace mymath {
	
	template<class T, size_t n>
	class vector {
		T values[n];
	public:
		using iterator = typename T*;
		using const_iterator = typename T const*;

		constexpr iterator begin() noexcept;
		constexpr iterator end() noexcept;

		constexpr const_iterator  begin() const noexcept;
		constexpr const_iterator end() const noexcept;

		vector(const T&);

		template<class U>
		vector(const vector<U, n>&);

		vector(const std::initializer_list<T>&);

		vector() = default;

		template<class U>
		vector<T, n>& copy(const vector<U, n>&);

		/* WARNING: Can change const object */
		static const vector<T, n>& fill(const vector<T, n>&, const T& elem = 0);

		template<class U>
		constexpr vector<T, n>& operator+=(const vector<U, n>&);

		template<class U>
		constexpr vector<T, n>& operator-=(const vector<U, n>&);

		template<class U>
		constexpr vector<T, n>& operator*=(const U&);

		template<class U>
		constexpr vector<T, n>& operator/=(const U&);
		
		constexpr T& operator[](size_t);

		constexpr T operator[](size_t) const;
	};

}

#include "../details/__macroses.inl"

namespace mymath {
	__GEN_MYMATH_VEC_PSEUDONYMS_09182(2, 2)
	__GEN_MYMATH_VEC_PSEUDONYMS_09182(3, 3)
	__GEN_MYMATH_VEC_PSEUDONYMS_09182(4, 4)
	__GEN_MYMATH_VEC_PSEUDONYMS_09182(5, 5)
}

