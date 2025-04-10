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


template<class BeginCond, class EndCond, class Func, class U0, class Ut0, class Uxx0, class StrClass>
struct conditions {
    BeginCond bc;
    EndCond ec;
    Func f;
    U0 u0;
    Ut0 ut0;
    Uxx0 uxx0;
	StrClass filename;
};

// each lambda has unique type so
template<class B, class E, class F, class U0, class Ut0, class Uxx0, class StrClass>
conditions<B, E, F, U0, Ut0, Uxx0, StrClass> make_conditions(B bc, E ec, F f, U0 u0, Ut0 ut0, Uxx0 uxx0, StrClass filename) {
    return conditions<B, E, F, U0, Ut0, Uxx0, StrClass>{ bc, ec, f, u0, ut0, uxx0, filename };
}

auto pr1 = make_conditions(
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                       				  					 // bc
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                 		    			  					 // ec
	[](double x, double t) {return 0; },         	  				// f
	[](double x) -> double {return std::sin(x * PI); },             // u0
	[](double x) { return 0; },             						// ut0
	[](double x) -> double {return -PI * PI * std::sin(x * PI); },  // uxx0
	"lab3_points1.txt"												// filename
);

auto pr2 = make_conditions(
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                       				  					 // bc
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                 		    			  					 // ec
	[](double x, double t) {return 0; },         	  				// f
	[](double x) -> double {return x * (x - 1); },          	    // u0
	[](double x) { return 0; },             						// ut0
	[](double x) -> double {return 2; },  							// uxx0
	"lab3_points_pr2.txt"											// filename
);

auto var10 = make_conditions(
	[](mymath::dynamic_vector<double>& line, double t = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },        				  					 // bc
	[](mymath::dynamic_vector<double>& line, double t) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0.5 *t; },     			  					 	   // ec
	[](double x, double t) {return 0; },	  						  // f
	[](double x) -> double {return (x + 1) * std::sin(x * PI); },     // u0
	[](double x) -> double {return x * (x + 1) ; },					  // ut0
	[](double x) -> double {return - PI * PI * (1 + x) * std::sin(PI * x) - 2 * PI * std::cos(PI * x); },  													   // uxx0
	"lab3_points_var10.txt"											  // filename
);

auto var22 = make_conditions(
	[](mymath::dynamic_vector<double>& line, double t = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0.5 *t; },										// bc
	[](mymath::dynamic_vector<double>& line, double t) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },													// ec
	[](double x, double t) {return 0; },								// f
	[](double x) -> double {return (2 - x) * std::sin(x * PI); },		// u0
	[](double x) -> double {return std::pow((x + 0.6), 2) ; },			// ut0
	[](double x) -> double {return -PI*PI * (2 - x) * std::sin(PI * x) - (PI + 1) * std::cos(PI * x); },  														// uxx0
	"lab3_points_var22.txt"												// filename
);

double f1(double x) {
	if ((x > -1./3.) && (x < 1./3.)) return 1;
	return 0;
}
double g2(double x) {
	if ((-1./2. < x) && (x < 1./2.)) return 1 - 2 * std::abs(x);
	return 0;
}

auto task1_pure = make_conditions(
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                       				  					 // bc
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                 		    			  					 // ec
	[](double x, double t) {return 0; },         	  				// f
	[](double x) -> double {return f1(x); },          	    		// u0
	[](double x) { return 0; },             						// ut0
	[](double x) -> double {return 0; },  							// uxx0
	"lab3_points_task1_pure.txt"									// filename
);

auto task2_pure = make_conditions(
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                       				  					 // bc
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                 		    			  					 // ec
	[](double x, double t) {return 0; },         	  				// f
	[](double x) -> double {return 0; },          	    			// u0
	[](double x) { return g2(x); },             					// ut0
	[](double x) -> double {return 0; },  							// uxx0
	"lab3_points_task2_pure.txt"									// filename
);

auto task3_pure = make_conditions(
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = std::sin(t); },                       				  					 // bc
	[](pddvec& line, double t =0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 0; },                 		    			  					 // ec
	[](double x, double t) {return 0; },         	  				// f
	[](double x) -> double {return 0; },          	    			// u0
	[](double x) { return 0; },             						// ut0
	[](double x) -> double {return 0; },  							// uxx0
	"lab3_points_task3_pure.txt"									// filename
);

template<class Type>
struct parameters {
	size_t n;
	Type L;
	Type a;
	Type step;
	Type tau;
	Type bt; // begin time
	Type start_point;
	std::string filename;
};

auto standard_param = parameters<double>{
	100,			// n
	2.,				// L
	1.,				// a
	2./(100-1),		// step
	2./(100-1)/10,	// tau
	2./(100-1)/10,	// bt
	0
};

auto pr2_param = parameters<double>{
	100,			// n
	1.,				// L
	1.,				// a
	1./(100-1),		// step
	1./(100-1)/10,	// tau
	1./(100-1)/10,	// bt
	0
};

auto task1_1 = parameters<double>{
	100,			// n
	4,				// L
	1.,				// a
	4./(100-1),		// step
	4./(100-1)/10,	// tau
	4./(100-1)/10,	// bt
	-2.,			// start point
	"lab3_points_task1_pure1.txt"
};

auto task1_2 = parameters<double>{
	100,			// n
	4,				// L
	5.,				// a
	4./(100-1),		// step
	4./(100-1)/10,	// tau
	4./(100-1)/10,	// bt
	-2.,			// start point
	"lab3_points_task1_pure2.txt"
};

auto task1_3 = parameters<double>{
	100,			// n
	4,				// L
	7.5,			// a
	4./(100-1),		// step
	4./(100-1)/10,	// tau
	4./(100-1)/10,	// bt
	-2.,			// start point
	"lab3_points_task1_pure3.txt"
};

auto task1_4 = parameters<double>{
	100,			// n
	4,				// L
	10.,			// a
	4./(100-1),		// step
	4./(100-1)/10,	// tau
	4./(100-1)/10,	// bt
	-2.,			// start point
	"lab3_points_task1_pure4.txt"
};

auto task2_1 = parameters<double>{
	100,			// n
	2,				// L
	1.,				// a
	2./(100-1),		// step
	2./(100-1)/10,	// tau
	2./(100-1)/10,	// bt
	-1.,			// start point
	"lab3_points_task2_pure1.txt"
};

auto task2_2 = parameters<double>{
	100,			// n
	2,				// L
	5.,				// a
	2./(100-1),		// step
	2./(100-1)/10,	// tau
	2./(100-1)/10,	// bt
	-1.,			// start point
	"lab3_points_task2_pure2.txt"
};

auto task2_3 = parameters<double>{
	100,			// n
	2,				// L
	7.5,			// a
	2./(100-1),		// step
	2./(100-1)/10,	// tau
	2./(100-1)/10,	// bt
	-1.,			// start point
	"lab3_points_task2_pure3.txt"
};

auto task2_4 = parameters<double>{
	100,			// n
	2,				// L
	10.,			// a
	2./(100-1),		// step
	2./(100-1)/10,	// tau
	2./(100-1)/10,	// bt
	-1.,			// start point
	"lab3_points_task2_pure4.txt"
};

auto task3_1 = parameters<double>{
	100,					// n
	4. * PI,				// L
	1.,						// a
	4. * PI/(100-1),		// step
	4. * PI/(100-1)/10,		// tau
	4. * PI/(100-1)/10,		// bt
	0.,						// start point
	"lab3_points_task3_pure1.txt"
};

auto task3_2 = parameters<double>{
	100,					// n
	4. * PI,				// L
	5.,						// a
	4. * PI/(100-1),		// step
	4. * PI/(100-1)/10,		// tau
	4. * PI/(100-1)/10,		// bt
	0.,						// start point
	"lab3_points_task3_pure2.txt"
};

auto task3_3 = parameters<double>{
	100,					// n
	4. * PI,				// L
	7.5,					// a
	4. * PI/(100-1),		// step
	4. * PI/(100-1)/10,		// tau
	4. * PI/(100-1)/10,		// bt
	0.,						// start point
	"lab3_points_task3_pure3.txt"
};

auto task3_4 = parameters<double>{
	100,					// n
	4. * PI,				// L
	10.,					// a
	4. * PI/(100-1),		// step
	4. * PI/(100-1)/10,		// tau
	4. * PI/(100-1)/10,		// bt
	0.,						// start point
	"lab3_points_task3_pure4.txt"
};

auto test_param = parameters<double>{
	100,				// n
	1./2.,				// L
	1.,					// a
	1./2./(100-1),		// step
	1./2./(100-1)/10,	// tau
	1./2./(100-1)/10,	// bt
	-2.
};


int main() {
	auto conds1 = task3_pure;
	auto params1 = task3_4;

	std::ofstream file1(conds1.filename, std::ios::trunc);
	file1.close(); // clear file
	
	mymath::wave_scheme<double, decltype(conds1.bc), decltype(conds1.ec), decltype(conds1.f)> ws(params1.bt, params1.tau, params1.n, params1.step, params1.a, conds1.bc, conds1.ec, conds1.f);
	
	ws.start_point = params1.start_point;
	ws.init(conds1.u0, conds1.ut0, conds1.uxx0);

	for (int k = 0; k < 5000; ++k) {
		ws.next();
		ws.points_out(params1.filename);
	}

	std::cout << "\nyo motherfucker\n";

	return 0;
}

