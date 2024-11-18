#include "mymath/mymath.h"
#include <list>
#include <iostream>
#include <vector>

using test_T = double;
constexpr const size_t SIZE = 4;
constexpr const double EPS = 1e-5;
using point = mymath::dvec2;
using pvec = mymath::dynamic_vector<point>;
using dvec = std::vector<double>;
constexpr const double PI = 3.1415926535;

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


void poly_fit_uniform(kaganvec& vec, size_t n, double l, double r) {
	auto labs = r - l;
	for (size_t i = 0; i != n - 1; ++i){
		
		add_point(vec, {l + i * labs / n, ex_f4(l + i * labs / n)});
	}
	add_point(vec, { r, ex_f4(r)});
}

void poly_fit_chebi(kaganvec& vec, size_t n, double l, double r) {
	auto labs = r - l;
	for (size_t i = 0; i != n + 1; ++i) {
		test_T x = (l + r) / 2 + ((r - l) / 2 ) * std::cos((2 * i + 1) * PI / (2*(n + 1)));
		add_point(vec, { x, ex_f4(x) });
	}

}

// spline related starts here
dvec progonka(const dvec& a, const dvec& b, const dvec& c, const dvec& d) 
{
    // Находим число n - количество неизветсных
    size_t n = b.size();

    // X - решение системы
    dvec X(n, 0);

    // Задаем коэффициенты прогонки
    dvec alpha(n, 0);
    dvec beta(n, 0);

    // Находим первые коэффициенты прогонки из 1-го уравнения
    alpha[1] = c[0] / b[0];
    beta[1] = d[0] / b[0];

	test_T denom;
    // Находим остальные коэффициенты прогонки
    for (size_t i = 1; i < n - 1; i++) {
        denom = b[i] - a[i] * alpha[i];
        alpha[i + 1] = c[i] / denom;
        beta[i + 1] = (d[i] + a[i] * beta[i]) / denom;
    }

    // Находим X(n) - последний элемент вектора решения
    X[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 1]) / (b[n - 1] - a[n - 1] * alpha[n - 1]);

    // Находим все остальные элементы вектора решения
    for (int i = n - 2; i >= 0; i--) {
        X[i] = alpha[i + 1] * X[i + 1] + beta[i + 1];
    }

    return X;
}

// Вычисляющая функция (вектор значений сплайна на заданной сетке)
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

        // Если поданные значения лежат вне сетки, то экстраполируем крайними полиномами сплайна:
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

// Построение полиномов (вычисление a, b, c, d)

dvec spline(const dvec& mesh, const dvec& F, const dvec& x) {
    size_t n = mesh.size();

    // Initialize arrays
    dvec h(n, 0);
    dvec g(n, 0);

    // Calculate h and g
    for (size_t i = 1; i < n; ++i) {
        h[i] = mesh[i] - mesh[i-1];
        g[i] = (F[i] - F[i-1]) / h[i];
    }

    // нижняя, главная, верхняя диагональ
    dvec A(n-2, 0);
    dvec B(n-2, 0);
    dvec C(n-2, 0);
	// правая часть
    dvec D(n-2, 0);

    // граничные условия на вторую производную
	// их значения можно менять, чтобы точность интерполяции была выше
    test_T diff2_a = 0;
    test_T diff2_b = 0;

	// члены будущего решения прогонкой
	// c0 = 0
    test_T c1 = diff2_a / 2;
    test_T cn_minus_2 = diff2_b / 2;
	// cn = 0 (именно n, размерность n+1, c[n-1])
	// cn_minus_2 по факту c[n-1]

    // Инициализация (A[0] равно A_2 в человеческом понимании, A[1] если бы начинали с 0.)
	// Не считаем A[1] по человечески, потому что сразу записываем решение c[0] = 0 и c[1], => эти коэффициенты не нужны.
	A[0] = 0;
    B[0] = -2 * (h[1] + h[2]);
    C[0] = h[2];
    D[0] = -(3 * (g[2] - g[1]) - c1 * h[1]);

	// Как бы A[2], A_3
    for (size_t i = 1; i < n-3; ++i) {
        B[i] = -2 * (h[i+1] + h[i+2]);
        A[i] = h[i+1];
        C[i] = h[i+2];
        D[i] = -3 * (g[i+2] - g[i+1]);
    }

	// Как бы B[n-2] и A[n-2], просто c[n-1] вычислили, тоже знаем решение поэтому не считаем
    B[n-3] = -2 * (h[n-2] + h[n-1]);
    A[n-3] = h[n-2];
	C[n-3] = 0;

	// для случая, когда размерность 3 коэффициенты в матрице прогонки примут другое значение
    if (n == 3) {
        D[n-3] = -(3 * (g[n] - g[n-1]) - cn_minus_2 * h[n] - c1 * h[1]);
    } else {
        D[n-3] = -(3 * (g[n] - g[n-1]) - cn_minus_2 * h[n]);
    }

    // подгон под нужную систему
	// на этой хуйне все ломается
	// for (size_t i = 0; i < n; i++) {
	// 	B[i] = -B[i];
	// 	D[i] = -D[i];
	// }

    dvec c_old = progonka(A, B, C, D);

    // дополнить c[0] = 0, c1 и c[n-1] ()
	dvec c(n, 0);
    c[1] = c1;
	for (size_t i = 2; i < n-1; i++) {
		c[i] = c_old[i-2];
	}
	c[n-1] = cn_minus_2;

    // Calculate remaining coefficients
    dvec a(n, 0);
    dvec b(n, 0);
    dvec d(n, 0);

    for (size_t i = 1; i < n-1; ++i) {
        a[i] = F[i-1];
        b[i] = g[i] - (c[i+1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i+1] - c[i]) / (3 * h[i]);
    }
	// для c_{n+1} = 0, i = n -1:
	a[n - 1] = F[n - 2];
	b[n - 1] = g[n-1] - 2 * c[n-1] * h[n - 1] / 3;
	d[n - 1] = - c[n - 1] / (3 * h[n - 1]);

    // Return the spline function
    return spline_eval(x, n, mesh, a, b, c, d);
}


int main(){
	/*
	kaganvec j;

	size_t n = 100;
	double l = -1, r = 1;
	double labs = r - l;
	poly_fit_chebi(j, 10, -1, 1);
	
	for (size_t i = 0; i != n + 1; ++i) {
		std::cout << '{' << l + i * labs / n << ", " << evaluate(j, l + i * labs / n) << "}, \n";
	}
	*/
	std::cout << "добрый день.\n";
	//dvec mesh = {-1, -0.5, 0, 0.5, 1};
	dvec mesh = {-1, -0.5, 0, 0.5, 1};
	dvec x = {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1};
	dvec F = {1, 2, 3, 4, 5};

	dvec sol = spline(mesh, F, x);
	size_t n = sol.size();
	for (size_t i = 0; i < n; i++) {
		std::cout << "{ " << sol[i] << ", " << x[i] << " }\n";
	}
	return 0;
}