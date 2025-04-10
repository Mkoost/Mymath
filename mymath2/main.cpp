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
	auto bc = [](mymath::dynamic_vector<double>& line, double t  = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; };
	auto f = [](double x, double t) {return 0; };


	auto u0 = [](double x) -> double {return std::sin(x * PI); };
	auto ut0 = [](double x) -> double {return 0; };
	auto uxx0 = [](double x) -> double {return - PI * PI * std::sin(x * PI); };
	
	auto u02 = [](double x) -> double {return x * (x - 1); };
	auto ut02 = [](double x) -> double {return 0; };
	auto uxx02 = [](double x) -> double {return 2; };
	
	auto u03 = [](double x) -> double {return (x + 1) * std::cos(x * PI); };
	auto ut03 = [](double x) -> double {return x * (x + 1) ; };
	auto uxx03 = [](double x) -> double {return - PI * PI * (1 + x) * std::cos(PI * x) - 2 * PI * std::sin(PI * x); };
	auto bc3 = [](mymath::dynamic_vector<double>& line, double t = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; };
	auto ec3 = [](mymath::dynamic_vector<double>& line, double t = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0.5 * t; };
	
	size_t n = 100;
	double L = 2;
	double bt = 0,  step = L / n,  a = 1., tau = step / (a * 10);
	

	mymath::wave_scheme<double, decltype(bc), decltype(bc), decltype(f)> ws(bt, tau, n, step, a, bc, bc, f);
	ws.init(u0, ut0, uxx0);
	ws.next(10);

	return 0;
}

