#pragma once
#include "../settings.h"
#include "../inline/dynamic_vector.inl.h"

namespace mymath {
	template<class T>
	struct difference_scheme {
		using __difference_scheme_K_coef_func = T(*)(T _tau, T _h);
		using __difference_scheme_bc_func = T(*)(T tau);

		using __difference_scheme_bc_approx = T(*)(
			T tau,
			dynamic_vector<T> elems,
			dynamic_vector<T> h_elems,
			__difference_scheme_bc_func bc, 
			__difference_scheme_K_coef_func K
			);

		using __difference_scheme_algo = T(*)(
			T _tau,
			T bt,
			size_t _h, // kolichestvo uzlov internal
			__difference_scheme_bc_func bbc, // begin boundary condition
			__difference_scheme_bc_func ebc, // end boundary condition
			__difference_scheme_bc_approx bbc_approx,
			__difference_scheme_bc_approx ebc_approx,
			dynamic_vector<T>* pl; // prev_layer
			dynamic_vector<T>* nl, // new_layer
			); 

		__difference_scheme_bc_func begin_bc = nullptr;
		__difference_scheme_bc_func end_bc = nullptr;
		
		T tau;
		size_t h;
		T end_time;

		dynamic_vector<T> next_layer;
		dynamic_vector<T> prev_layer;
		

		/*bool next();
		void log(std::string out = "");
		void out_layer(std::string out = "");
		
		void gen_grid(T h_, T tau_) {
			
		};
<<<<<<< HEAD
		void set_begin_time(T begin_time);*/
=======
>>>>>>> 7d5dff1ea53a6a8e55eb0afb3cd2c7b73bd5ebc7
	};

}