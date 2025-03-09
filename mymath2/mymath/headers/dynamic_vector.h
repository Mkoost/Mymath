#pragma once
#include<initializer_list>

namespace mymath {
	template<class T>
	class dynamic_vector {
		T* values = nullptr;
		size_t n_ = 0;

	public:
		// TODO
		using iterator = T*;
		using const_iterator = T const*;

		dynamic_vector() = default;

		// TODO
		dynamic_vector(const T&, size_t n);

		dynamic_vector(const std::initializer_list<T>&);

		// TODO
		template<class U>
		dynamic_vector(const dynamic_vector<U>& A);

		dynamic_vector(const dynamic_vector<T>& A);

		dynamic_vector(dynamic_vector<T>&& A);

		~dynamic_vector();

		// TODO
		template<class U>
		dynamic_vector<T>& copy(const dynamic_vector<U>&);

		dynamic_vector<T>& move(dynamic_vector<T>&) ;

		dynamic_vector<T>& move(T*, size_t);

		size_t columns() const;
		size_t rows() const;
		size_t size() const;
		unsigned id() const;

		// TODO
#if 0
		// TODO
		constexpr iterator begin()  noexcept;
		constexpr iterator end() noexcept;

		// TODO
		constexpr const_iterator  begin() const noexcept;
		constexpr const_iterator end() const noexcept;
#endif
		// TODO
#if 0
		reverse_iterator rbegin();
		reverse_iterator rend();
#endif



		// TODO
#if 0
		reverse_const_iterator rbegin() const;
		reverse_const_iterator rend() const;
#endif

		// TODO
		/* WARNING: Can change const object */
		static const dynamic_vector<T>& fill(const dynamic_vector<T>&, const T& elem = 0);

		// TODO
		constexpr dynamic_vector<T>& operator=(const dynamic_vector<T>&);

		// TODO
		constexpr dynamic_vector<T>& operator=(dynamic_vector<T>&&);

		// TODO
		constexpr dynamic_vector<T>& operator+=(const dynamic_vector<T>&);

		// TODO
		constexpr dynamic_vector<T>& operator-=(const dynamic_vector<T>&);

		// TODO
		constexpr dynamic_vector<T>& operator*=(const T&);

		// TODO
		constexpr dynamic_vector<T>& operator/=(const T&);

		// TODO
		constexpr T& operator[](size_t);

		// TODO
		constexpr T operator[](size_t) const;

	};

}