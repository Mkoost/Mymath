#pragma once

#include "../matrix.inl.h"
#include "../vector.inl.h"
#include "../dynamic_matrix.inl.h"
#include "../dynamic_vector.inl.h"
#include "../../details/__expr.inl"
#include "../../settings.h"
#include <cmath>
#include <list>

namespace mymath {
	
	template<class T, class U>
	std::list<mymath::dynamic_vector<T>> runge_kutta_4(T bt, T et, const dynamic_vector<T>& init, const dynamic_vector<U>& f, T step, T stepmin = 1e-13, double eps=1e-8, size_t iter=4294967296) {
		if (init.size() != f.size())
			throw(std::invalid_argument("Init vector size and f size are different"));
		if (bt > et)
			throw(std::invalid_argument("Begin time is greater then end time"));

		const size_t n = init.size();

		std::list<dynamic_vector<T>> res;
		dynamic_vector<T> yn(init);
		res.push_back(yn);
		
		dynamic_vector<T> k(0, n);
		dynamic_vector<T> ki(0, n);
		dynamic_vector<T> tmp(0, n);
		do {
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
			res.push_back(yn);
			if (et - bt < step) step = et - bt;
		} while (bt <= et && step > stepmin);
		std::cout << bt << '\n';
		return res;

	}



}