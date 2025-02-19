#pragma once

#include "../matrix.inl.h"
#include "../vector.inl.h"
#include "../dynamic_matrix.inl.h"
#include "../dynamic_vector.inl.h"
#include "../../details/__expr.inl"
#include "../../settings.h"
#include <cmath>
#include <list>
#include <utility>


namespace mymath {
	

	namespace{
		template<class T, class U>
		struct __special_de_func_for_explict_euler {
			size_t state = 0;
			T tau = 0;
			U* f = nullptr;
			dynamic_vector<T>* yi = nullptr;

			__special_de_func_for_explict_euler<T, U>& operator[](size_t st) {
				state = st;
				return *this;
			}
	

			T operator() (const T& t, const mymath::dynamic_vector<T>& y) {
				return y[state] - (*yi)[state] - tau * (*f)[state](t, y);
			}

			size_t size() const { return yi->size(); };
		};

		template<class T, class U>
		struct __special_de_func_for_symmetrical_scheme {
			size_t state = 0;
			T tau = 0;
			U* f = nullptr;
			dynamic_vector<T>* yi = nullptr;
			dynamic_vector<T>* fi = nullptr;
			__special_de_func_for_symmetrical_scheme<T, U>& operator[](size_t st) {
				state = st;
				return *this;
			}


			T operator() (const T& t, const mymath::dynamic_vector<T>& y) {
				return 2*(y[state] - (*yi)[state]) - tau * ((*f)[state](t, y) + (*fi)[state]);
			}

			size_t size() const { return yi->size(); };
		};

	}

	template<class T, class U>
	std::list<mymath::dynamic_vector<T>> runge_kutta_4(T bt, T et, const dynamic_vector<T>& init, const U& f, T step, T stepmin = 1e-13, double eps=1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		const size_t n = init.size();

		std::list<dynamic_vector<T>> res;
		dynamic_vector<T> yn(init);
		res.push_back(yn);
		size_t iter = 0;
		dynamic_vector<T> k(0, n);
		dynamic_vector<T> ki(0, n);
		dynamic_vector<T> tmp(0, n);
		do {
			++iter;
			for (size_t i = 0; i < n; ++i)
				ki[i] = f[i](bt, yn);
			
			k += ki;
			tmp = yn + (step / 2) * ki;
			bt += step / 2;

			for (size_t i = 0; i < n; ++i)
				ki[i] = f[i](bt, tmp);
			
			tmp = yn + (step / 2) * ki;
			ki *= 2;
			k += ki;

			for (size_t i = 0; i < n; ++i)
				ki[i] = f[i](bt, tmp);

			tmp = yn + step * ki;

			ki *= 2;
			k += ki;
			bt += step / 2;

			for (size_t i = 0; i < n; ++i)
				ki[i] = f[i](bt, tmp);

			k += ki;
			k /= 6;

			yn += step * k;
			if (iter % save_iter == 0) {
				res.push_back(yn);
			}
			if (et - bt < step) step = et - bt;
		} while ((bt <= et) && (step > stepmin) && (iter != max_iter));
		return res;

	}


	template<class T, class U>
	std::list<mymath::dynamic_vector<T>> runge_kutta_2(T bt, T et, const dynamic_vector<T>& init, const U& f, T step, T stepmin = 1e-13, double eps = 1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		const size_t n = init.size();

		std::list<dynamic_vector<T>> res;
		dynamic_vector<T> yn(init);
		res.push_back(yn);
		size_t iter = 0;
		dynamic_vector<T> k(0, n);
		dynamic_vector<T> ki(0, n);
		dynamic_vector<T> tmp(0, n);
		do {
			++iter;
			for (size_t i = 0; i < n; ++i)
				ki[i] = f[i](bt, yn);
			k += ki;
			
			tmp = yn + step * ki;
			for (size_t i = 0; i < n; ++i)
				ki[i] = f[i](bt, tmp);
			k += ki;

			bt += step;			

			k /= 2;

			yn += step * k;

			if (iter % save_iter == 0)
				res.push_back(yn);
			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		std::cout << bt << '\n';
		return res;

	}
	
	template<class T, class U>
	std::list<std::pair<T, dynamic_vector<T>>> runge_kutta_4_autostep_fast(T bt, T et, const dynamic_vector<T>& init, const U& f, T step, T stepmin = 1e-13, double eps = 1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		const size_t n = init.size();

		std::list<std::pair<T, dynamic_vector<T>>> res;
		dynamic_vector<T> yn(init);

		std::pair<T, dynamic_vector<T>> elem;
		elem.first = bt;
		elem.second = yn;
		res.push_back(elem);

		size_t iter = 0;
		dynamic_vector<T> k(0, n);
		dynamic_vector<T> ki(0, n);
		dynamic_vector<T> tmp(0, n);



		do {
			++iter;
			double eps_actual;
			tmp = yn;
			auto bigstep = runge_kutta_4(bt, et, tmp, f, step, stepmin, eps, 1, 1);
			auto smallstep = runge_kutta_4(bt, et, tmp, f, step / 2, stepmin, eps, 2, 2);
			auto bg_b = bigstep.end();
			auto bg_s = smallstep.end();
			if (bigstep.size() < 2 || bigstep.size() < 2) break;
			--bg_b;
			--bg_s;
			eps_actual = cube_norm((*bg_b) - (*bg_s)) / 15;
			step = std::pow(eps / eps_actual, 0.2) * step;

			smallstep = runge_kutta_4(bt, et, tmp, f, step / 2, stepmin, eps, 2, 2);
			bg_s = smallstep.begin();
			++bg_s;
			yn = (*bg_s);

			bt += step;

			if (iter % save_iter == 0) {
				elem.first = bt;
				elem.second = yn;
				res.push_back(elem);
				std::cout << bt << '\n';
			}

			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		return res;

	}

	template<class T, class U>
	std::list<std::pair<T, dynamic_vector<T>>> runge_kutta_4_autostep_slow(T bt, T et, const dynamic_vector<T>& init, const U& f, T step, T stepmin = 1e-13, double eps = 1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));
		
		const size_t n = init.size();

		std::list<std::pair<T, dynamic_vector<T>>> res;
		dynamic_vector<T> yn(init);

		std::pair<T, dynamic_vector<T>> elem;
		elem.first = bt;
		elem.second = yn;
		res.push_back(elem);

		size_t iter = 0;
		dynamic_vector<T> k(0, n);
		dynamic_vector<T> ki(0, n);
		dynamic_vector<T> tmp(0, n);

		do {
			++iter;
			double eps_actual;
			tmp = yn;
			auto bigstep = runge_kutta_4(bt, et, tmp, f, step, stepmin, eps, 1, 1);
			auto smallstep = runge_kutta_4(bt, et, tmp, f, step / 2, stepmin, eps, 2, 2);
			auto bg_b = bigstep.end();
			auto bg_s = smallstep.end();
			if (bigstep.size() < 2 || bigstep.size() < 2) break;
			--bg_b;
			--bg_s;
			size_t iner_iter = 0;
			do {
				++iner_iter;
				eps_actual = cube_norm((*bg_b) - (*bg_s)) / 15;

				if (std::abs(eps_actual / eps) < 10 && std::abs(eps_actual / eps) > 0.1) { yn = (*bg_b);  break; }
				step *= std::pow(eps / eps_actual, 0.2);

				*bg_b = (*bg_s);
				
				smallstep = runge_kutta_4(bt, et, tmp, f, step / 2, stepmin, eps, 2, 2);
				bg_s = smallstep.end();
				--bg_s;
				if (smallstep.size() < 2) { yn = (*bg_b); break; }
			} while (true);
			
			if (smallstep.size() == 0 && iner_iter < 5) break;

			bt += step;

			if (iter % save_iter == 0) {
				elem.first = bt;
				elem.second = yn;
				res.push_back(elem);
				std::cout << bt << '\n';
			}
			
			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		return res;

	}


	template<class T, class U>
	std::list<mymath::dynamic_vector<T>> explicit_euler(T bt, T et, const dynamic_vector<T>& init, const U& f, T step, T stepmin = 1e-13, double eps = 1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		const size_t n = init.size();

		std::list<dynamic_vector<T>> res;
		dynamic_vector<T> yn(init);
		res.push_back(yn);
		size_t iter = 0;
		dynamic_vector<T> k(0, n);
		dynamic_vector<T> tmp(0, n);
		do {
			++iter;
			for (size_t i = 0; i < n; ++i)
				k[i] = f[i](bt, yn);

			bt += step;
			yn += step * k;

			if (iter % save_iter == 0)
				res.push_back(yn);
			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		std::cout << bt << '\n';
		return res;

	}

	template<class T, class U>
	std::list<mymath::dynamic_vector<T>> implict_euler(T bt, T et, const dynamic_vector<T>& init, U& f, T step, T stepmin = 1e-13, double eps = 1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		const size_t n = init.size();

		std::list<dynamic_vector<T>> res;
		dynamic_vector<T> yn(init);
		res.push_back(yn);
		size_t iter = 0;
		__special_de_func_for_explict_euler<T, U> method;
		method.tau = step;
		method.f = &f;
		method.yi = &yn;
		do {
			++iter;
			yn = eq_system_solve(bt, yn, method, eps);
			bt += step;

			if (iter % save_iter == 0)
				res.push_back(yn);
			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		return res;

	}

	template<class T, class U>
	std::list<mymath::dynamic_vector<T>> symmetrical_scheme(T bt, T et, const dynamic_vector<T>& init, U& f, T ststep, T stepmin = 1e-13, double eps = 1e-8, size_t save_iter = 1, size_t max_iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		T step = ststep;
		const size_t n = init.size();

		size_t iter = 0;
		std::list<dynamic_vector<T>> res;
		dynamic_vector<T> yn(init);
		dynamic_vector<T> fn(0, yn.size());
		res.push_back(yn);

		__special_de_func_for_symmetrical_scheme<T, U> method;
		method.tau = step;
		method.f = &f;
		method.yi = &yn;
		method.fi = &fn;

		do {
			++iter;
			for (size_t i = 0; i < yn.size(); ++i)
				fn[i] = f[i](bt, yn);
			yn = eq_system_solve(bt, yn, method, eps);
			bt += step;

			if (iter % save_iter == 0)
				res.push_back(yn); 

			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		return res;

	}

	/*template<class T, class U>
	std::list<mymath::dynamic_vector<T>> predictor_corrector(T bt, T et, const dynamic_vector<T>& init, U& f, T ststep, T stepmin = 1e-13, double eps = 1e-8, size_t iter = 4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		T step = ststep;
		dynamic_vector<T> yn_corr(init);
		dynamic_vector<T> yn_pred(init);
		dynamic_vector<T> fn_pred(init);
		std::list<dynamic_vector<T>> res = runge_kutta_4(bt, bt + ststep / 10, init, f, ststep / 40);
		

		auto bg = res.begin();
		auto fi = bg;
		end += 3;

		T pred_step = ststep / 10;
		do {
			dynamic_vector<T>::fill(yn_pred, 0);
			fi = bg;

			yn_pred += 55 * (*fi);
			++fi;
			yn_pred -= 59 * (*fi);
			++fi;
			yn_pred += 37 * (*fi);
			++fi;
			yn_pred -= 9 * (*fi);
			
			yn_pred /= 24;
			yn_pred *= pred_step;
			yn_pred += yn_corr;
			++bg;
			
			fn_pred()

			fi = bg;

			bt += step;
			pred_step = step;
			res.push_back(yn);
			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		return res;

	}*/
}