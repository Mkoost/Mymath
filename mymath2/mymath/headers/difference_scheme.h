#pragma once
#include "../settings.h"
#include "../inline/dynamic_vector.inl.h"
#include "../inline/dynamic_matrix.inl.h"
#include "../inline/vector.inl.h"
#include "../inline/utilities.inl"
namespace mymath {
	namespace { 
		constexpr const double __mixed_difference_scheme_sigma = 0.;
	
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
		T step;
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
				mymath::utilities::print(prev_layer);
				begin_time += tau;
			}
		}


		static void mixed_algo(difference_scheme* ds) {
			dynamic_vector<T> line(0, 4);
			T sigma = __mixed_difference_scheme_sigma;

			// left cond

			ds->bbc_approx(ds, line);

			for (size_t i = 0; i < 3; ++i)
				ds->progonka_buff[0][i] = line[i];

			ds->next_layer[0] = line[3];

			// right cond


			ds->ebc_approx(ds, line);

			for (size_t i = 0; i < 3; ++i)
				ds->progonka_buff[ds->progonka_buff.rows() - 1][i] = line[i];

			ds->next_layer[ds->next_layer.size() - 1] = line[3];

			T h = 1.0 / (ds->prev_layer.size() - 1);
			T tau = ds->tau;

			for (size_t i = 1; i < ds->prev_layer.size() - 1; ++i) {
				T a1 = ds->K(0, h * (i - 0.5));
				T a2 = ds->K(0, h * (i + 0.5));
				T w1 = a1 * (ds->prev_layer[i] - ds->prev_layer[i - 1]) / h;
				T w2 = a2 * (ds->prev_layer[i + 1] - ds->prev_layer[i]) / h;
				ds->progonka_buff[i][0] = sigma / h * a1;
				ds->progonka_buff[i][1] = -(sigma / h * (a1 + a2) + ds->cp_coef * h / tau);
				ds->progonka_buff[i][2] = sigma / h * a2;
				ds->next_layer[i] = -(ds->cp_coef * h / tau * ds->prev_layer[i] + (1 - sigma) * (w2 - w1));

			}

			diag3_solver(ds->progonka_buff, ds->next_layer);

		};
		
		static void semi_explicit_algo(difference_scheme* ds) {
			dynamic_vector<T> line(0, 4);
			
			// left cond

			ds->bbc_approx(ds, line);

			for (size_t i = 0; i < 3; ++i)
				ds->progonka_buff[0][i] = line[i];

			ds->next_layer[0] = line[3];

			// right cond


			ds->ebc_approx(ds, line);

			for (size_t i = 0; i < 3; ++i)
				ds->progonka_buff[ds->progonka_buff.rows() - 1][i] = line[i];

			ds->next_layer[ds->next_layer.size() - 1] = line[3];

			T h = ds->step;
			T tau = ds->tau;

			for (size_t i = 1; i < ds->prev_layer.size() - 1; ++i) {
				T dk1 = (ds->K((ds->prev_layer[i] + ds->prev_layer[i - 1]) / 2 + 1e-10, 0) - ds->K((ds->prev_layer[i] + ds->prev_layer[i - 1]) / 2, 0)) / 1e-10;
				T dk2 = (ds->K((ds->prev_layer[i + 1] + ds->prev_layer[i]) / 2 + 1e-10, 0) - ds->K((ds->prev_layer[i + 1] + ds->prev_layer[i]) / 2, 0)) / 1e-10;
				T a1 = ds->K((ds->prev_layer[i] + ds->prev_layer[i - 1]) / 2, 0) + tau * dk1;
				T a2 = ds->K((ds->prev_layer[i + 1] + ds->prev_layer[i]) / 2, 0) + tau * dk2;

				/*T a1 = 0.5 * (ds->K(ds->prev_layer[i - 1], h * (i - 1)) + ds->K(ds->prev_layer[i], h * i));
				T a2 = 0.5 * (ds->K(ds->prev_layer[i], h * i) + ds->K(ds->prev_layer[i + 1], h * (i + 1)));*/


				ds->progonka_buff[i][0] = -tau * a1;
				ds->progonka_buff[i][1] = h * h + tau * (a1 + a2);
				ds->progonka_buff[i][2] = -tau * a2;
				ds->next_layer[i] = h * h * ds->prev_layer[i];

			}

			/*utilities::print(ds->progonka_buff);
			utilities::print(ds->next_layer);*/
			diag3_solver(ds->progonka_buff, ds->next_layer);

		};

		static void semi_explicit_algo_simple_iter(difference_scheme* ds) {
			dynamic_vector<T> prev_layer(ds->prev_layer);
			dynamic_vector<T> line(0, 4);

			double norm = 0;
			size_t iter = 0;
			do{
				++iter;

				// left cond

				ds->bbc_approx(ds, line);

				for (size_t i = 0; i < 3; ++i)
					ds->progonka_buff[0][i] = line[i];

				ds->next_layer[0] = line[3];

				// right cond


				ds->ebc_approx(ds, line);

				for (size_t i = 0; i < 3; ++i)
					ds->progonka_buff[ds->progonka_buff.rows() - 1][i] = line[i];

				ds->next_layer[ds->next_layer.size() - 1] = line[3];

				T h = ds->step;
				T tau = ds->tau;

				for (size_t i = 1; i < ds->prev_layer.size() - 1; ++i) {


					T a1 = 0.5 * (ds->K(prev_layer[i - 1], h * (i - 1)) + ds->K(ds->prev_layer[i], h * i));
					T a2 = 0.5 * (ds->K(prev_layer[i], h * i) + ds->K(prev_layer[i + 1], h * (i + 1)));


					ds->progonka_buff[i][0] = -tau * a1;
					ds->progonka_buff[i][1] = h * h + tau * (a1 + a2);
					ds->progonka_buff[i][2] = -tau * a2;
					ds->next_layer[i] = h * h * ds->prev_layer[i];

				}

				/*utilities::print(ds->progonka_buff);
				utilities::print(ds->next_layer);*/
				diag3_solver(ds->progonka_buff, ds->next_layer);
				
				norm = 0;
				for (size_t i = 0; i < ds->next_layer.size(); ++i)
					norm = std::max(norm, std::abs(ds->next_layer[i] - prev_layer[i]));
				
				
			} while (norm > 1e-10);
			std::cout << iter << "\n";
		};

		dynamic_vector<T> next_layer;
		dynamic_vector<T> prev_layer;
		dynamic_matrix<T> progonka_buff;
		
	};

	
	

	template<class Type_, class BeginCond_, class EndCond_, class F_>
	class wave_scheme {
	public:
		dynamic_vector<Type_> next_layer;
		dynamic_vector<Type_> prev_layer;
		dynamic_vector<Type_> pprev_layer;
		dynamic_matrix<Type_> progonka_buff;

		Type_ tau;
		size_t k_tau = 1;
		Type_ begin_time;
		Type_ step;
		size_t n;
		Type_ start_point;

		Type_ a;
		F_ f;
		BeginCond_ bc;
		EndCond_ ec;



		wave_scheme(Type_ begin_time_, Type_ tau_, size_t n_, Type_ step_, Type_ a_, BeginCond_ bc_, EndCond_ ec_, F_ f_) : 
			begin_time(begin_time_),
			tau(tau_),
			n(n_),
			step(step_),
			a(a_),
			bc(bc_),
			ec(ec_),
			f(f),
			next_layer(0, n_),
			prev_layer(0, n_),
			pprev_layer(0, n_),
			progonka_buff(0, n_, 3),
			start_point(0.0){};

		

		static void cross_scheme(wave_scheme<Type_, BeginCond_, EndCond_, F_>& ws){
			dynamic_vector<Type_> line(0, 4);
			
			// left cond
			auto n = ws.n;

			ws.bc(line, ws.begin_time);

			for (size_t i = 0; i < 3; ++i)
				ws.progonka_buff[0][i] = line[i];

			ws.next_layer[0] = line[3];

			ws.ec(line, ws.begin_time);

			for (size_t i = 0; i < 3; ++i)
				ws.progonka_buff[n - 1][i] = line[i];

			ws.next_layer[n - 1] = line[3];

			Type_ h = ws.step;
			Type_ tau = ws.tau;

			for (size_t i = 1; i < n - 1; ++i) {
				ws.progonka_buff[i][1] = 1.;
				ws.next_layer[i] = 2 * ws.prev_layer[i] - ws.pprev_layer[i] + std::pow((ws.a * tau) / h, 2) * (ws.prev_layer[i + 1] - 2 * ws.prev_layer[i] + ws.prev_layer[i - 1]) + ws.f(ws.start_point + i * h, ws.begin_time);
			}

			diag3_solver(ws.progonka_buff, ws.next_layer);
		}

		void next(size_t k = 1) {
			k_tau += k;
			std::cout << progonka_buff.columns() << " " << progonka_buff.rows() << "\n";
			for (size_t i = 0; i < k; ++i) {
				cross_scheme(*this);
				dynamic_vector<Type_> tmp;
				tmp.move(pprev_layer);
				pprev_layer.move(prev_layer);
				prev_layer.move(next_layer);
				next_layer.move(tmp);
				begin_time += tau;
			}
		};

		template<class f1, class f2, class f3>
		void init(f1 u, f2 ut, f3 uxx) {
			for (int i = 0; i < n; ++i) pprev_layer[i] = u(start_point + step * i);
			for (int i = 0; i < n; ++i) prev_layer[i] = pprev_layer[i] + tau * ut(start_point + step * i) + uxx(start_point + step * i) * a * a * tau * tau / 2;
		}

		void points_out(std::string filename) {
			std::ofstream file1(filename, std::ios::app);
			for (size_t i = 0; i < n; ++i) {
				file1 << k_tau * tau << " " << start_point + i*step << " " << next_layer[i] << "\n";
			}

			file1.close();
		}
	};

}


//
//
//
//
