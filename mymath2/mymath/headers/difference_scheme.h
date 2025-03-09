#pragma once
#include "../settings.h"
#include "../inline/dynamic_vector.inl.h"
#include "../inline/dynamic_matrix.inl.h"
#include "../inline/vector.inl.h"
#include "../inline/utilities.inl"
namespace mymath {
	namespace { constexpr const double __mixed_difference_scheme_sigma = 1;
	
	template<class T>
	void diag3_solver(dynamic_matrix<T>& buff, dynamic_vector<T> & X)
	{
		size_t n = X.size();


		dynamic_vector<T> alpha(0.0, n);
		dynamic_vector<T> beta(0.0, n);

		alpha[1] = -buff[0][2] / buff[0][1];
		beta[1] = X[0] / buff[0][1];

		T denom;
		for (size_t i = 1; i < n - 1; i++) {
			denom = buff[i][1] + buff[i][0] * alpha[i];
			alpha[i + 1] = -buff[i][2] / denom;
			beta[i + 1] = (X[i] - buff[i][0] * beta[i]) / denom;
		}

		X[n - 1] = (X[n - 1] - buff[n - 1][0] * beta[n - 1]) / (buff[n - 1][1] + buff[n - 1][0] * alpha[n - 1]);

		for (int i = n - 2; i >= 0; i--) {
			X[i] = alpha[i + 1] * X[i + 1] + beta[i + 1];
		}

	}
	
	
	}




	template<class T>
	class difference_scheme{
	public:
		using difference_scheme_K_coef_func = T(*)(T _u, T _x);
		using difference_scheme_bc_func = T(*)(T tau);
		using difference_scheme_bc_approx = void(*)(
			difference_scheme* ds,
			dynamic_vector<T>& pl
			);

		using difference_scheme_algo = void(*)(difference_scheme* ds);

		T tau;
		T begin_time;
		T cp_coef;
		difference_scheme_bc_approx bbc_approx;
		difference_scheme_bc_approx ebc_approx;
		difference_scheme_K_coef_func K;
		difference_scheme_algo algo;

		difference_scheme(
			T _tau, 
			T _start_time,
			size_t _h,
			dynamic_vector<T>& start_u0,
			T _cp_coef, 
			difference_scheme_bc_approx _bbc, 
			difference_scheme_bc_approx _ebc, 
			difference_scheme_K_coef_func _K,
			difference_scheme_algo _algo){

			_h += 2;
			tau = _tau;
			begin_time = _start_time;
			cp_coef = _cp_coef;
			bbc_approx = _bbc;
			ebc_approx = _ebc;
			K = _K;
			algo = _algo;

			next_layer.move(new T[_h], _h);
			prev_layer.copy(start_u0);
			progonka_buff.move(new T[_h * 3], _h, 3);
		}


		void next(size_t k = 1) {
			for (auto i = 0; i < k; ++i) {
				algo(this);
				dynamic_vector<T> tmp;
				tmp.move(prev_layer);
				prev_layer.move(next_layer);
				next_layer.move(tmp);
				
				begin_time += tau;
			}
		}


		static void mixed_algo(difference_scheme* ds) {
			dynamic_vector<T> line(0, 4);
			T sigma = __mixed_difference_scheme_sigma;

			// left cond
			for (size_t i = 0; i < 2; ++i)
				line[i] = ds->prev_layer[i];

			ds->bbc_approx(ds, line);

			for (size_t i = 0; i < 3; ++i)
				ds->progonka_buff[0][i] = line[i];

			ds->next_layer[0] = line[3];

			// right cond
			for (size_t i = 0; i < 2; ++i)
				line[i] = ds->prev_layer[ds->prev_layer.size() - 2 + i];

			ds->ebc_approx(ds, line);

			for (size_t i = 0; i < 3; ++i)
				ds->progonka_buff[ds->progonka_buff.rows() - 1][i] = line[i];

			ds->next_layer[ds->next_layer.size() - 1] = line[3];

			T h = 1.0 / (ds->prev_layer.size());
			T tau = ds->tau;
			
			for (size_t i = 1; i < ds->prev_layer.size() - 1; ++i) {
				T a1 = 1 / ds->K(0, h * (i - 0.5));
				T a2 = 1 / ds->K(0, h * (i + 0.5));
				T w1 = a1 * (ds->prev_layer[i] - ds->prev_layer[i - 1]) / h;
				T w2 = a2 * (ds->prev_layer[i + 1] - ds->prev_layer[i]) / h;
				ds->progonka_buff[i][0] = sigma / h * a1;
				ds->progonka_buff[i][1] = -(sigma / h * (a1 + a2) + ds->cp_coef * h / tau);
				ds->progonka_buff[i][2] = sigma / h * a2;
				ds->next_layer[i] = -(ds->cp_coef * h / tau * ds->prev_layer[i] + (1 - sigma) * (w2 - w1));

			}

			diag3_solver(ds->progonka_buff, ds->next_layer);

		}
		
		dynamic_vector<T> next_layer;
		dynamic_vector<T> prev_layer;
		dynamic_matrix<T> progonka_buff;
		
		

	};

	


}


//
//
//
//
