#include "mymath/mymath.h"
#include <list>
#include <iostream>
#include <vector>
#include <fstream>
#include <functional>

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

// разделить bc и ec на horizontal и vert
template<class BeginCond_hor, class EndCond_hor, class BeginCond_vert, class EndCond_vert, class Func, class U0, class StrClass>
struct conditions_2D {
    BeginCond_hor bc_h;
    EndCond_hor ec_h;
	BeginCond_vert bc_v;
    EndCond_vert ec_v;
    Func f;
    U0 u0;
	StrClass filename;
};

// each lambda has unique type so
template<class Bh, class Eh, class Bv, class Ev, class F, class U0, class StrClass>
conditions_2D<Bh, Eh, Bv, Ev, F, U0, StrClass> make_conditions(Bh bc_h, Eh ec_h, Bv bc_v, Ev ec_v, F f, U0 u0, StrClass filename) {
    return conditions_2D<Bh, Eh, Bv, Ev, F, U0, StrClass>{ bc_h, ec_h, bc_v, ec_v, f, u0, filename };
}

// (edit line[3] to set boundary conditions in both cases)
auto pr1 = make_conditions(
	[](pddvec& line, double t = 0, double ui = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                       				  					 // bc_h
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                 		    			  					 // ec_h
	
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                       				  					 // bc_v
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                 		    			  					// ec_v
	
	[](double x, double y) {return 0; },         	  		// f
	[](double x, double y) -> double {return 0; },             		// u0
	"lab4_points1.txt"												// filename
);

// здесь будем передавать названия файла
// будет генериться три файла с таким названием:
// <name>_t.txt, <name>_x.txt, <name>yt.txt, 
// step_x = L / (n-1), step_x = L / (n-1)
// tau = min(step_{x, y}) / 10, bt = tau
template<class Type>
struct parameters {
	size_t n;
	size_t m;
	Type L;
	Type step_x;
	Type step_y;
	Type tau;
	Type bt; // begin time
	Type start_point_x;
	Type start_point_y;
	std::string filename;
};

auto standard_param = parameters<double>{
	100,			// n
	100,			// m
	1.,				// L
	1./(100-1),		// step_x
	1./(100-1),		// step_y
	1./(100-1)/10,	// tau
	1./(100-1)/10,	// bt
	0,				// start_point_x
	0,				// start_point_y
	"lab4_pr1"
};

/////////////////////////////////////////////
int main() {
	auto conds1 = pr1;
	auto params1 = standard_param;

	std::ofstream file1(params1.filename + "_t.txt", std::ios::trunc);
	std::ofstream file2(params1.filename + "_points.txt", std::ios::trunc);
	file1.close();
	file2.close(); // clear file
	
	mymath::lapl2d_scheme<double, 
						  decltype(conds1.bc_h), 
						  decltype(conds1.ec_h),
		                  decltype(conds1.bc_v),
		                  decltype(conds1.ec_v),
		                  decltype(conds1.f)> 
				ws(params1.bt,
					params1.tau,
					params1.n,
					params1.m,
					params1.step_x,
					params1.step_y,
					conds1.bc_h,
					conds1.ec_h,
					conds1.bc_v,
					conds1.ec_v,
					conds1.f);
	
	ws.st_x = params1.start_point_x;
	ws.st_y = params1.start_point_y;
	ws.init(conds1.u0);

	for (int k = 0; k < 1; ++k) {
		ws.next();
		ws.state_out(params1.filename);
	}


	std::cout << "\nyo motherfucker\n";

	return 0;
}

