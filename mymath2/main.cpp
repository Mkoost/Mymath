#include "mymath/mymath.h"
#include <list>
#include <iostream>

using test_T = double;
constexpr const size_t SIZE = 4;
constexpr const double EPS = 1e-5;
using point = mymath::dvec2;
using pvec = mymath::dynamic_vector<point>;
constexpr const double PI = 3.1415926535;

struct kaganvec {
	std::list<test_T> fs;
	std::list<point> ps;
};

void add_point(kaganvec& vec, point p) {
	if (vec.ps.size()) {
		vec.ps.push_back(p);
		test_T f = 0;
		auto psb = vec.ps.begin();
		auto psi = vec.ps.begin();
		++psi;
		auto scnd = psi;
		auto pse = vec.ps.end();
		for (size_t i = 0; i != vec.ps.size()-1; ++i) {
			test_T tmp = 1;
			++psb;
			while (psi != pse) {
				if (psi != psb)
					tmp *= (*psb)[0] - (*psi)[0];
				++psi;
			}
			f += (*psb)[1] / tmp;
			psi = scnd;
		}
		f -= (*vec.fs.rbegin());
		f /= p[0] - (*vec.ps.begin())[0];
		vec.fs.push_back(f);
	}
	else {
		vec.fs.push_back(p[1]);
		vec.ps.push_back(p);
	}
};

test_T evaluate(kaganvec& vec, test_T x) {
	auto ps = vec.ps.begin();
	
	test_T res = (*ps)[1];
	test_T tmp = 1;
	auto f = vec.fs.begin();
	++f;
	for (auto endf = vec.fs.end(); f != endf; ++f) {
		tmp *= (x - (*ps)[0]);
		res += tmp * (*f);
		++ps;
	}
	return res;

}


test_T ex_f1(test_T x) {
	return  x * x;
}

test_T ex_f2(test_T x) {
	return 1 / (1 + x * x);
}

test_T ex_f3(test_T x) {
	return 1 / std::atan(1 + 10 * x * x);
}

test_T ex_f4(test_T x) {
	return std::pow((4 * x * x * x + 2 * x * x - 4 * x + 2), std::sqrt(2))
		+ std::sin(1/(5 + x - x * x)) - 5;
}

test_T ex_f5(test_T x) {
	return std::pow(2.71828182846, x);
}



void poly_fit_uniform(kaganvec& vec, size_t n, double l, double r) {
	auto labs = r - l;
	for (size_t i = 0; i != n - 1; ++i){
		
		add_point(vec, {l + i * labs / n, ex_f5(l + i * labs / n)});
	}
	add_point(vec, { r, ex_f5(r)});
}

void poly_fit_chebi(kaganvec& vec, size_t n, double l, double r) {
	auto labs = r - l;
	for (size_t i = 0; i != n + 1; ++i) {
		test_T x = (l + r) / 2 + ((r - l) / 2 ) * std::cos((2 * i + 1) * PI / (2*(n + 1)));
		add_point(vec, { x, ex_f4(x) });
	}

}



int main(){
	kaganvec j;

	size_t n = 100;
	double l = -1, r = 1;
	double labs = r - l;
	poly_fit_uniform(j, 11, 0, 2);
	std::cout << 9.0250134994341209264717 - evaluate(j, 2.2);
	

	return 0;
}