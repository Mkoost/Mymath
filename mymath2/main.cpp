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
using eqSys = mymath::dynamic_vector<double(*)(const double, const mymath::dynamic_vector<double>)>;
struct kaganvec {
	std::list<test_T> fs;	// functions?
	std::list<point> ps;	// points?
};


/*
test_T func1(test_T x, test_T y) {
	return  x * x - y * y - 15;
}
test_T func2(test_T x, test_T y) {
	return  x * y + 4;
}
*/
/*
test_T func1(test_T x, test_T y) {
	return sqrt( sqrt(x)+sqrt(y) ) - 8 * (sqrt(x)-sqrt(y));
}
test_T func2(test_T x, test_T y) {
	return sqrt(x)+sqrt(y)-4;
}
*/


// spline related starts here
dvec progonka(const dvec& a, const dvec& b, const dvec& c, const dvec& d)
{
	// ������� ����� n - ���������� �����������
	size_t n = b.size();

	// X - ������� �������
	dvec X(n, 0);

	// ������ ������������ ��������
	dvec alpha(n, 0);
	dvec beta(n, 0);

	// ������� ������ ������������ �������� �� 1-�� ���������
	alpha[1] = c[0] / b[0];
	beta[1] = d[0] / b[0];

	test_T denom;
	// ������� ��������� ������������ ��������
	for (size_t i = 1; i < n - 1; i++) {
		denom = b[i] - a[i] * alpha[i];
		alpha[i + 1] = c[i] / denom;
		beta[i + 1] = (d[i] + a[i] * beta[i]) / denom;
	}

	// ������� X(n) - ��������� ������� ������� �������
	X[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 1]) / (b[n - 1] - a[n - 1] * alpha[n - 1]);

	// ������� ��� ��������� �������� ������� �������
	for (int i = n - 2; i >= 0; i--) {
		X[i] = alpha[i + 1] * X[i + 1] + beta[i + 1];
	}

	return X;
}

// ����������� ������� (������ �������� ������� �� �������� �����)
dvec spline_eval(const dvec& x, size_t n, const dvec& mesh,
	const dvec& a, const dvec& b,
	const dvec& c, const dvec& d) {
	// Pre-calculate the length of x
	size_t N = x.size();

	// Create an array for spline values
	dvec spl(N, 0);

	// Find spline values for each x(j)
	for (size_t j = 0; j < N; ++j) {

		size_t I = 0;
		// Find the interval in which x(j) is located
		for (size_t i = 0; i < n - 1; ++i) {
			if (x[j] >= mesh[i] && x[j] <= mesh[i + 1]) {
				I = i + 1;
				break;
			}
		}

		// ���� �������� �������� ����� ��� �����, �� �������������� �������� ���������� �������:
		if (x[j] < mesh[0]) {
			I = 1;
		}
		if (x[j] > mesh[n - 1]) {
			I = n;
		}

		// Calculate the spline value at x(j)
		spl[j] = a[I] + b[I] * (x[j] - mesh[I - 1]) + c[I] * pow(x[j] - mesh[I - 1], 2) + d[I] * pow(x[j] - mesh[I - 1], 3);
	}

	return spl;
}

// ���������� ��������� (���������� a, b, c, d)

dvec spline(const dvec& mesh, const dvec& F, const dvec& x) {
	size_t n = mesh.size();

	// Initialize arrays
	dvec h(n, 0);
	dvec g(n, 0);

	// Calculate h and g
	for (size_t i = 1; i < n; ++i) {
		h[i] = mesh[i] - mesh[i - 1];
		g[i] = (F[i] - F[i - 1]) / h[i];
	}

	// ������, �������, ������� ���������
	dvec A(n - 2, 0);
	dvec B(n - 2, 0);
	dvec C(n - 2, 0);
	// ������ �����
	dvec D(n - 2, 0);

	// ��������� ������� �� ������ �����������
	// �� �������� ����� ������, ����� �������� ������������ ���� ����
	test_T diff2_a = 0;
	test_T diff2_b = 0;

	// ����� �������� ������� ���������
	// c0 = 0
	test_T c1 = diff2_a / 2;
	test_T cn_minus_2 = diff2_b / 2;
	// cn = 0 (������ n, ����������� n+1, c[n-1])
	// cn_minus_2 �� ����� c[n-1]

	// ������������� (A[0] ����� A_2 � ������������ ���������, A[1] ���� �� �������� � 0.)
	// �� ������� A[1] �� �����������, ������ ��� ����� ���������� ������� c[0] = 0 � c[1], => ��� ������������ �� �����.
	A[0] = 0;
	B[0] = -2 * (h[1] + h[2]);
	C[0] = h[2];
	D[0] = -(3 * (g[2] - g[1]) - c1 * h[1]);

	// ��� �� A[2], A_3
	for (size_t i = 1; i < n - 3; ++i) {
		B[i] = -2 * (h[i + 1] + h[i + 2]);
		A[i] = h[i + 1];
		C[i] = h[i + 2];
		D[i] = -3 * (g[i + 2] - g[i + 1]);
	}

	// ��� �� B[n-2] � A[n-2], ������ c[n-1] ���������, ���� ����� ������� ������� �� �������
	B[n - 3] = -2 * (h[n - 2] + h[n - 1]);
	A[n - 3] = h[n - 2];
	C[n - 3] = 0;

	// ��� ������, ����� ����������� 3 ������������ � ������� �������� ������ ������ ��������
	if (n == 3) {
		D[n - 3] = -(3 * (g[n-1] - g[n - 2]) - cn_minus_2 * h[n-1] - c1 * h[1]);
	}
	else {
		D[n - 3] = -(3 * (g[n-1] - g[n - 2]) - cn_minus_2 * h[n-1]);
	}

	// ������ ��� ������ �������
	// �� ���� ����� ��� ��������
	// for (size_t i = 0; i < n; i++) {
	// 	B[i] = -B[i];
	// 	D[i] = -D[i];
	// }

	dvec c_old = progonka(A, B, C, D);

	// ��������� c[0] = 0, c1 � c[n-1] ()
	dvec c(n, 0);
	c[1] = c1;
	for (size_t i = 2; i < n - 1; i++) {
		c[i] = c_old[i - 2];
	}
	c[n - 1] = cn_minus_2;

	// Calculate remaining coefficients
	dvec a(n, 0);
	dvec b(n, 0);
	dvec d(n, 0);

	for (size_t i = 1; i < n - 1; ++i) {
		a[i] = F[i - 1];
		b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}
	// ��� c_{n+1} = 0, i = n -1:
	a[n - 1] = F[n - 2];
	b[n - 1] = g[n - 1] - 2 * c[n - 1] * h[n - 1] / 3;
	d[n - 1] = -c[n - 1] / (3 * h[n - 1]);

	// Return the spline function
	return spline_eval(x, n, mesh, a, b, c, d);
}




std::list<mymath::dvec2> roots_location(double a, double b, size_t n, double (*f)(double), double eps = EPS) {
	double step = (b - a) / n;
	std::list<mymath::dvec2> res;

	mymath::dvec2 v{0, 0};
	double f1 = f(a);
	
	if (f1 < 0)
		v = { a, 0 };
	else
		v = { 0, a };

	double x1 = a;

	for (size_t i = 1; i != n + 1; ++i) {
		
		double x2 = a + step * i;
		double f2 = f(x2);

		if (f1 == 0){
			if (f2 < 0){
				v[0] = x2;
				v[1] = x1;
				res.push_back(v);
			}
			else{
				v[1] = x2;
				v[0] = x1;
				res.push_back(v);
			}
		}
		else if (f1 * f2 > 0) {
			if (f1 < 0)
				v[0] = x2;
			else
				v[1] = x2;
			
		}
		else {
			if (f1 < 0){
				v[1] = x2;
				res.push_back(v);
			}
			else{
				v[0] = x2;
				res.push_back(v);
			}
		}

		f1 = f2;
		x1 = x2;

	}

	return res;

};

// tau -- шаг разностной сетки, с начального приближения в левом конце считаем все точки на области
std::list<pddvec> explicit_Euler(pddvec grid, pddvec y0, pddvec (*func)(double tn, pddvec yn)) {
	std::list<pddvec> res;
	res.push_back(y0);

	size_t n = grid.size();
	for (int i = 1; i < n; ++i) {
		res.push_back(res.back() + (grid[i] - grid[i - 1]) * func(grid[i], res.back()));
	}

	return res;
}
double Pr(double pr) {
	return 1;
}
double Pl(double pr) {
	return -1;
}

int main() {

	mymath::difference_scheme<double>::difference_scheme_bc_approx
		bbc_1type = [](mymath::difference_scheme<double>* ds, mymath::dynamic_vector<double>& vec) -> void {double u0 = std::pow(2. * 5. * 5. / 0.5, 0.5);  for (int i = 0; i < 4; ++i) vec[i] = 0; vec[1] = 1; vec[3] = u0 * std::pow(ds->begin_time, 0.5); },
		ebc_1type = [](mymath::difference_scheme<double>* ds, mymath::dynamic_vector<double>& vec) -> void {for (int i = 0; i < 4; ++i) vec[i] = 0; vec[1] = 1; vec[3] = 0; },
		ebc_2type = [](mymath::difference_scheme<double>* ds, mymath::dynamic_vector<double>& vec) -> void {
			for (int i = 0; i < 4; ++i) vec[i] = 0; 
			 
			double h  = 1. / ds->next_layer.size();
			double an = 1. / ds->K(0, h * (ds->next_layer.size() - 0.5));
			double sig_an_h = mymath::__mixed_difference_scheme_sigma * an / h;
			double cp_h_2tau = ds->cp_coef * h / (2 * ds->tau);
			double w = an * (ds->prev_layer[ds->prev_layer.size() - 1] - ds->prev_layer[ds->prev_layer.size() - 1]) / h;;
			double kappa = sig_an_h / (cp_h_2tau + sig_an_h);
			double mu = (cp_h_2tau * ds->prev_layer[ds->prev_layer.size() - 1] +
						mymath::__mixed_difference_scheme_sigma * Pr(ds->begin_time + ds->tau)
						+ (1 - mymath::__mixed_difference_scheme_sigma) * (Pr(ds->begin_time) - w)) / (cp_h_2tau + sig_an_h);
			vec[1] = 1;
			vec[0] = -kappa;
			vec[3] = mu; 
		},
		bbc_2type = [](mymath::difference_scheme<double>* ds, mymath::dynamic_vector<double>& vec) -> void {
		for (int i = 0; i < 4; ++i) vec[i] = 0;

		double h = 1. / ds->next_layer.size();
		double an = 1. / ds->K(0, h * 0.5);
		double sig_an_h = mymath::__mixed_difference_scheme_sigma * an / h;
		double cp_h_2tau = ds->cp_coef * h / (2 * ds->tau);
		double w = an * (ds->prev_layer[1] - ds->prev_layer[0]) / h;;
		double kappa = sig_an_h / (cp_h_2tau + sig_an_h);
		double mu = (cp_h_2tau * ds->prev_layer[0] +
			mymath::__mixed_difference_scheme_sigma * Pl(ds->begin_time + ds->tau)
			+ (1 - mymath::__mixed_difference_scheme_sigma) * (Pl(ds->begin_time) - w)) / (cp_h_2tau + sig_an_h);
		vec[1] = 1;
		vec[2] = -kappa;
		vec[3] = mu;
		};
	mymath::difference_scheme<double>::difference_scheme_K_coef_func 
		K = [](double u, double x) -> double { return 0.5 * std::pow(u,2); };

	size_t n = 798;
	mymath::dynamic_vector<double> u0(0, n + 2);
	double ttau = 0.0001;
	mymath::difference_scheme<double> ds(ttau, ttau, n, u0, 1.0, bbc_1type, ebc_1type, K, mymath::difference_scheme<double>::semi_explicit_algo);
	ds.step = 0.0125;
	size_t k = 15000;
	ds.next(k);

	for (size_t i = 0; i < n; ++i) u0[i] = std::sqrt(2. * 5. * 2. * (5. * k * ttau - i * ds.step));
	mymath::utilities::print(ds.prev_layer);
	std::cout << "\n";
	mymath::utilities::print(u0);
	std::cout << "\n";
	mymath::utilities::print(ds.prev_layer - u0);
	return 0; 
}
