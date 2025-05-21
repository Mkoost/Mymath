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
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0, double f_ij = 0, double x = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                       				  					 // bc_h
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0, double f_ij = 0, double x = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                 		    			  					 // ec_h
	
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0, double f_ij = 0, double x = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                       				  					 // bc_v
	[](pddvec& line, double t = 0, double ui_min_1 = 0, double u_i = 0, double ui_pl_1 = 0, double h1 = 0, double h2 = 0, double f_ij = 0, double x = 0) { line[0] = 0; line[1] = 1; line[2] = 0; line[3] = 1; },                 		    			  					// ec_v
	
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
/***************************************************/
// ЛАБА 5 ---------------------------------------

// ЧЕТНОЕ КОЛИЧЕСТВО УЗЛОВ n!!!!
template<class Ker_>
void fredgolm_quadrature_matrix(double a, double b, size_t n, double lambda, Ker_ K, mymath::dynamic_matrix<double>& A) {
	size_t N = n + 1;
	double h = (b - a) / n;
	mymath::dynamic_matrix<double> mat(0, N, N);
	mymath::dynamic_matrix<double>::diag(mat);

	for (size_t i = 0; i < N; ++i){

		mat[i][0] -= lambda * h / 2.0 * K(a + i * h, a );
		mat[i][N - 1] -= lambda * h / 2.0 * K(a + i * h, a + (N - 1) * h);
		for (size_t j = 1; j + 1 < N ; j += 1) { // 1/3 4/3 1/3
			mat[i][j] -= lambda * h  * K(a + i * h, a + j * h);
		}
		
	}

	A.move(mat);
}


// ЧЕТНОЕ КОЛИЧЕСТВО УЗЛОВ n!!!!
template<class Ker_, class Func_>
void  fredgolm_quadrature_solve(double a,
								double b,
								size_t n,
								double lambda,
								Ker_ K,
								Func_ f,
								mymath::dynamic_vector<double>& u) {
	size_t N = n + 1;
	double h = (b - a) / n;
	
	mymath::dynamic_matrix<double> mat;
	mymath::dynamic_vector<double> fs(0, n+1);
	fredgolm_quadrature_matrix(a, b, n, lambda, K, mat);
	

	for (size_t i = 0; i < N; ++i)
		fs[i] = f(a + i*h);
	mymath::data_structs::base_data_dynamic_vector_matrix<double, 1, 1> res = mymath::qr_solve(mat, fs);
	u.move(res.vec[0]);
	/*auto res = mymath::relax_iteration(mat, fs, 1e-8);
	u.move(res);*/
}


template<class Ker_, class Func_>
void fredgolm_simple_iter_estim(mymath::dynamic_vector<double>& u, 
								  double a, double b, size_t n,
								  double lambda,
								  Ker_ K,
								  Func_ f,
								  mymath::dynamic_vector<double>& ans,
								  double eps = 1e-6) {
	

	double h = (b - a) / double(n);
	int N = n + 1;
	double norm = 0;
	mymath::dynamic_vector<double> fx(0, N);

	for (size_t i = 0; i < N; ++i)
		fx[i] = f(a + i * h);

	mymath::dynamic_vector<double> u_res = u;
	mymath::dynamic_vector<double> u_new = fx;
	 do{
		 norm = 0;
		for (int i = 0; i < N; ++i){
			double tmp = 0;
			for (int j = 0; j + 1 < N; ++j)
				tmp += h * (u_res[j] * K(a + i * h, a + j * h) + u_res[j + 1] * K(a + i * h, a + (j + 1) * h)) / 2.0;
			u_new[i] += lambda * tmp;
			norm = std::max(norm,std::abs(u_new[i] - u_res[i]));
		}
		std::swap(u_new, u_res);
		u_new = fx;
	 } while (norm > eps);
	ans.move(u_res);
}


template<class Func_>
void contour_solve(
	size_t n,
	Func_ f,
	mymath::dynamic_vector<double>& u) {

	auto Qx = [](double x, double y, double xi, double eta) { return -(y - eta) / ((x - xi) * (x - xi) + (y - eta) * (y - eta)) / (2.0 * PI); };
	auto Qy = [](double x, double y, double xi, double eta) { return (x - xi) / ((x - xi) * (x - xi) + (y - eta) * (y - eta)) / (2.0 * PI); };
	mymath::dynamic_matrix<double> mat(0, n + 1, n + 1);
	mymath::dynamic_vector<double> fs(0, n + 1);
	double lj = 2.0 * PI / n;

	for (int i = 0; i < n; ++i) {
		
		double kxi = std::cos(2.0 * PI * (i - 0.5) / n);
		double kyi = std::sin(2.0 * PI * (i - 0.5) / n);

		double nxi = kxi;
		double nyi = kyi;

		double lj = 2.0 * PI / n;
		
		fs[i] = f(kxi, kyi);

		for (int j = 0; j < n; ++j) {
			double cxj = std::cos(2.0 * PI * (i - 1.0) / n);
			double cyj = std::sin(2.0 * PI * (i - 1.0) / n);

			mat[i][j] = (nxi * Qx(kxi, kyi, cxj, cxj) + nyi * Qy(kxi, kyi, cxj, cxj)) * lj;
		}
		mat[i][n] = 1.0;
	}

	for (int i = 0; i < n; ++i) mat[n][i] = lj;
	
	mymath::data_structs::base_data_dynamic_vector_matrix<double, 1, 1> res = mymath::qr_solve(mat, fs);
	u.move(res.vec[0]);
};



int main() {
	int n = 1000;
	mymath::dynamic_vector<double> uinit(0, n + 1);
	mymath::dynamic_vector<double> u;

	auto K1 = [](double x, double s) {return 1.0; };
	auto f1 = [](double x) {return 1.0; };
	auto lambda1 = 0.25;
	auto u1 = [](double x) -> double {return 2.; };

	auto K2 = [](double x, double s) {return 1 - x * std::cos(x * s); };
	auto f2 = [](double x) {return (1.0 + std::sin(x)) / 2.0; };
	auto lambda2 = 0.5;

	auto K3 = [](double x, double s) {return 1 - x * std::cos(x * s); };
	auto f3 = [](double x) {return x * x + std::sqrt(x); };
	auto lambda3 = 0.5;

	auto K4 = [](double x, double s) {return std::exp(s + x); };
	auto f4 = [](double x) {return x; };
	auto res4 = [](double x) {return (-2 * std::exp(x) - 3 * x + std::exp(2)*x) / (-3 + std::exp(2)); };
	auto lambda4 = 1.0;
	double a4 = 0.0, b4 = 1.0;

	double a2 = 0, b2 = 1.0;
	double a3 = 0.1, b3 = 2.0;

	fredgolm_quadrature_solve(a2, b2, n, lambda2, K2, f2, u); // 2
	
	for (int i = 0; i < n + 1; ++i)
		uinit[i] = res4(a4 + i * (b4 - a4) / n);

	double err = std::abs(1.0 - u[0]);
	for (int i = 1; i < n + 1; ++i) err = std::max(err, std::abs(1.0 - u[i]));
	
	std::cout << err;


	// Пример 1 из методички
	//5: 0.0014739
	//10: 0.000367309
	//20: 9.17545e-05
	//50: 1.46775e-05
	//100: 3.66925e-06
	//1000: 3.66921e-08
	return 0;
}

