#include "mymath/mymath.h"
#include <list>
#include <iostream>
#include <vector>
#include <fstream>

using test_T = double;
constexpr const size_t SIZE = 4;
constexpr const double EPS = 1e-6;
constexpr const double D_EPS = 1e-8;
using point = mymath::dvec2;
using pvec = mymath::dynamic_vector<point>;
using dvec = std::vector<double>;
using pddvec = mymath::dynamic_vector<double>;
constexpr const double PI = 3.1415926535897932;
constexpr const double E = 2.71828182845904523536;
using eqSys = mymath::dynamic_vector<double(*)(const double, const mymath::dynamic_vector<double>)>;
struct kaganvec {
	std::list<test_T> fs;	// functions?
	std::list<point> ps;	// points?
};



// 
int main() {
	auto bc = [](mymath::dynamic_vector<double>& line) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; };
	auto f = [](double x, double t) {return 0; };
	size_t n = 100;
	double bt = 0, tau = 0.1, step = PI / n,  a = 1.;
	

	mymath::wave_scheme ws(bt, tau, n, step, a, bc, bc, f);
	for (int i = 0; i < n; ++i) ws.pprev_layer[i] = std::sin(step * i * PI);

	return 0;
}

