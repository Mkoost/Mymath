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
constexpr const double PI = 3.1415926535897932;
using eqSys = mymath::dynamic_vector<double(*)(double, mymath::dynamic_vector<double>)>;
struct kaganvec {
	std::list<test_T> fs;	// functions?
	std::list<point> ps;	// points?
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
		for (size_t i = 0; i != vec.ps.size() - 1; ++i) {
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
	return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

test_T ex_f2(test_T x) {
	return std::sqrt(x + 1) - 1;
}

test_T ex_f3(test_T x) {
	return  35 * x * x * x - 67 * x * x - 3 * x + 3;
}

test_T ex_f4(test_T x) {
	return std::pow((4 * x * x * x + 2 * x * x - 4 * x + 2), std::sqrt(2))
		+ std::sin(1 / (5 + x - x * x)) - 5;
}

test_T ex_f5(test_T x) {
	return (x - 1) * (x - 1) * (x - 1);
}
test_T ex_f6(test_T x) {
	return x * x - 1;
}

test_T func(test_T x) {
	return  ex_f3(x);
}
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
test_T func1(test_T x, test_T y) {
	return pow(pow(x, 2) + 2 * pow(y, 2), 2) - 7 * (pow(x, 2) - 2 * pow(y, 2));
}
test_T func2(test_T x, test_T y) {
	return 3 * pow(x, 2) + 4 * pow(y, 2) - 14;
}
void poly_fit_uniform(kaganvec& vec, size_t n, double l, double r) {
	auto labs = r - l;
	for (size_t i = 0; i != n - 1; ++i) {
		add_point(vec, { l + i * labs / n, func(l + i * labs / n) });
	}
	add_point(vec, { r, func(r) });
}

void poly_fit_chebi(kaganvec& vec, size_t n, double l, double r) {
	for (size_t i = 1; i != n + 1; ++i) {
		test_T x = (l + r) / 2 + ((r - l) / 2) * std::cos((2 * i - 1) * PI / (2 * (n)));
		add_point(vec, { x, func(x) });
	}

}

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

test_T err_norm(kaganvec& j, size_t n, double l, double r) {
	test_T mx = 0;
	test_T labs = r - l;
	for (size_t i = 0; i < n + 1; ++i) {
		mx = std::max(mx, std::fabs(func(l + i * labs / n) - evaluate(j, l + i * labs / n)));
	}
	return mx;
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



/*
int main() {
	eqSys a(nullptr, 2);
	a[0] = [](double t, mymath::dynamic_vector<double> x) -> double {return 2 * x[0] + x[1] * x[1] - 1; };
	a[1] = [](double t, mymath::dynamic_vector<double> x) -> double {return 6 * x[0] - x[1] * x[1] + 1; };
	mymath::dynamic_vector<double> init = { 0.0, 0.0 };
	auto res = mymath::runge_kutta_4(0.0, 1.0, init, a, 1e-3);
	std::cout << res.size();
	return 0; 
}
*/

int main() {
	eqSys fl(nullptr, 2);
	fl[0] = [](double t, mymath::dynamic_vector<double> x) -> double {return x[0] * x[0] - x[1] * x[1] - 15; };
	fl[1] = [](double t, mymath::dynamic_vector<double> x) -> double {return x[0] * x[1] + 4; };
	eqSys fr(nullptr, 2);
	fr[0] = [](double t, mymath::dynamic_vector<double> x) -> double {return 0; };
	fr[1] = fr[0];

	mymath::dynamic_vector<double> st_p = { -3, 1.5 };


	auto x =  mymath::eq_system_solve(0.0, st_p, fl, fr);

	std::cout << x[0] << " " << x[1] << '\n';
}