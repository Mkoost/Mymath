#pragma once
#include <initializer_list>
namespace mymath {

	template<class T>
	class dynamic_matrix{

		// TODO
		class row {
			T* ptr;
		public:
			row(T*);
			T& operator[](size_t);
		};

		// TODO
		class const_row {
			T const* ptr;
		public:
			const_row(T const *);
			T operator[](size_t) const;
		};

		T* values = nullptr;
		size_t n_ = 0, m_ = 0;

	public:
		// TODO
		using iterator =  T*;
		using const_iterator =  T const*;

		dynamic_matrix() = default;

		// TODO
		dynamic_matrix(const T&, size_t n, size_t m);

		dynamic_matrix(const std::initializer_list<std::initializer_list<T>>&);

		// TODO
		template<class U>
		dynamic_matrix(const dynamic_matrix<U>& A);

		dynamic_matrix(const dynamic_matrix<T>& A);

		dynamic_matrix(dynamic_matrix<T>&& A);

		~dynamic_matrix();

		// TODO
		template<class U>
		dynamic_matrix<T>& copy(const dynamic_matrix<U>&);

		dynamic_matrix<T>& move(dynamic_matrix<T>&) ;
		dynamic_matrix<T>& move(T*, size_t n, size_t m) ;

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
		static const dynamic_matrix<T>& fill(const dynamic_matrix<T>&, const T& elem = 0);

		// TODO
		/* WARNING: Can change const object */
		static const dynamic_matrix<T>& diag(const dynamic_matrix<T>&, const T& elem = 1);

		// TODO
		constexpr dynamic_matrix<T>& operator=(const dynamic_matrix<T>&);
		
		// TODO
		constexpr dynamic_matrix<T>& operator=(dynamic_matrix<T>&&);

		// TODO
		constexpr dynamic_matrix<T>& operator+=(const dynamic_matrix<T>&);

		// TODO
		constexpr dynamic_matrix<T>& operator-=(const dynamic_matrix<T>&);

		// TODO
		constexpr dynamic_matrix<T>& operator*=(const T&);

		// TODO
		constexpr dynamic_matrix<T>& operator/=(const T&);

		// TODO
		constexpr row operator[](size_t);

		// TODO
		constexpr const_row operator[](size_t) const;


	};


}