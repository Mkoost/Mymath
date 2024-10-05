#include "mymath/mymath.h"
#include <iostream>


using test_T = double;
constexpr const size_t SIZE = 4;
constexpr const double EPS = 1e-10;

template<class T, class tmp_T>
using solver_ptr = mymath::data_structs::base_data_dynamic_vector_matrix<T, 1, 1>(*)(mymath::dynamic_matrix<T>&, mymath::dynamic_vector<T>&, tmp_T);

template<class T, class tmp_T>
using iter_solver_ptr = mymath::dynamic_vector<T>(*)(const mymath::dynamic_vector<T>&,
	const mymath::dynamic_matrix<T>&,
	tmp_T, mymath::dynamic_vector<T>&);


template<class T, class tmp_T>
mymath::data_structs::base_data_dynamic_vector_matrix<T, 1, 1> gauss(mymath::dynamic_matrix<T>& A, mymath::dynamic_vector<T>& b, tmp_T zero) {
	return mymath::gauss_solve(A, b, zero);
};

template<class T, class tmp_T>
mymath::data_structs::base_data_dynamic_vector_matrix<T, 1, 1> qr(mymath::dynamic_matrix<T>& A, mymath::dynamic_vector<T>& b, tmp_T zero) {
	return mymath::qr_solve(A, b, zero);
};

template<class T, class tmp_T>
mymath::dynamic_matrix<T>& jacobi(const mymath::dynamic_vector<T>&A,
	const mymath::dynamic_matrix<T>&b,
	tmp_T eps, mymath::dynamic_vector<T>& start) {
	return mymath::jacobi_iteration(A, b, eps, start);
};


template<class T, class tmp_T>
void test_solve(const std::string& mat_path, const std::string& vec_path, solver_ptr<T, tmp_T> solver) {
	mymath::dynamic_matrix<T> A = mymath::utilities::input_matrix_NxN<T>(mat_path);

	mymath::dynamic_matrix<T> A_inv(0, A.rows(), A.columns());

	mymath::dynamic_vector<T> b = mymath::utilities::input_vector_N<T>(vec_path);

	mymath::dynamic_matrix<T> B;

	mymath::dynamic_vector<T> x;



	try {
		auto some = solver(A, b, 1e-15);

		B.move(some.mat[0]);
		x.move(some.vec[0]);

	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what() << "\n\n";
	}
	catch (std::bad_alloc& e) {
		std::cerr << "Caught bad_alloc: " << e.what() << "\n\n";
	}

	std::cout << "Matrix A: \n";
	mymath::utilities::print(A);
	std::cout << "Matrix B: \n";
	mymath::utilities::print(B);
	std::cout << "Vector b: \n";
	mymath::utilities::print(b);
	std::cout << "\n";

	if (x.size() != 0) {

		auto y = mymath::multiply(A, x);
		std::cout << "Vector x: \n";
		mymath::utilities::print(x);
		std::cout << "\n";

		std::cout << "Vector y = Ax: \n";
		mymath::utilities::print(y);
		std::cout << "\n";

		std::cout << "|| b - y || = " << mymath::norm(b - y);
		std::cout << "\n\n\n";

		mymath::dynamic_vector<T> line(0, A.rows());

		line[0] = 1;
		for (size_t i = 0; i != line.size() - 1; ++i) {
			auto tmp = solver(A, line, 1e-15);

			for (size_t j = 0; j != line.size(); ++j)
				A_inv[j][i] = tmp.vec[0][j];

			std::swap(line[i], line[i + 1]);
		}


		auto tmp = solver(A, line, 1e-15);
		for (size_t j = 0; j != line.size(); ++j)
			A_inv[j][line.size() - 1] = tmp.vec[0][j];

		std::cout << "Inverse A: \n";
		mymath::utilities::print(A_inv);
		std::cout << "\n";

		std::cout << "|| A ||_1 = " << mymath::oct_norm(A) << "\n";
		std::cout << "|| A ^ -1 ||_1 = " << mymath::oct_norm(A_inv) << "\n";
		std::cout << "cond_1 A = " << mymath::oct_norm(A) * mymath::oct_norm(A_inv) << "\n";

		std::cout << "\n";
		
		std::cout << "|| A ||_inf = " << mymath::cube_norm(A) << "\n";
		std::cout << "|| A ^ -1 ||_inf = " << mymath::cube_norm(A_inv) << "\n";
		std::cout << "cond_inf A = " << mymath::cube_norm(A) * mymath::cube_norm(A_inv) << "\n";


	}
	else std::cout << "nope\n";


}

template<class T, class tmp_T>
void test_iter_solve(const std::string& mat_path, const std::string& vec_path, solver_ptr<T, tmp_T> solver) {
	mymath::dynamic_matrix<T> A = mymath::utilities::input_matrix_NxN<T>(mat_path);

	mymath::dynamic_matrix<T> A_inv(0, A.rows(), A.columns());

	mymath::dynamic_vector<T> b = mymath::utilities::input_vector_N<T>(vec_path);

	mymath::dynamic_matrix<T> B;

	mymath::dynamic_vector<T> x;



	try {
		auto some = solver(A, b, 1e-15);

		B.move(some.mat[0]);
		x.move(some.vec[0]);

	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what() << "\n\n";
	}
	catch (std::bad_alloc& e) {
		std::cerr << "Caught bad_alloc: " << e.what() << "\n\n";
	}

	std::cout << "Matrix A: \n";
	mymath::utilities::print(A);
	std::cout << "Matrix B: \n";
	mymath::utilities::print(B);
	std::cout << "Vector b: \n";
	mymath::utilities::print(b);
	std::cout << "\n";

	if (x.size() != 0) {

		auto y = mymath::multiply(A, x);
		std::cout << "Vector x: \n";
		mymath::utilities::print(x);
		std::cout << "\n";

		std::cout << "Vector y = Ax: \n";
		mymath::utilities::print(y);
		std::cout << "\n";

		std::cout << "|| b - y || = " << mymath::norm(b - y);
		std::cout << "\n\n\n";

		mymath::dynamic_vector<T> line(0, A.rows());

		line[0] = 1;
		for (size_t i = 0; i != line.size() - 1; ++i) {
			auto tmp = solver(A, line, 1e-15);

			for (size_t j = 0; j != line.size(); ++j)
				A_inv[j][i] = tmp.vec[0][j];

			std::swap(line[i], line[i + 1]);
		}


		auto tmp = solver(A, line, 1e-15);
		for (size_t j = 0; j != line.size(); ++j)
			A_inv[j][line.size() - 1] = tmp.vec[0][j];

		std::cout << "Inverse A: \n";
		mymath::utilities::print(A_inv);
		std::cout << "\n";

		std::cout << "|| A ||_1 = " << mymath::oct_norm(A) << "\n";
		std::cout << "|| A ^ -1 ||_1 = " << mymath::oct_norm(A_inv) << "\n";
		std::cout << "cond_1 A = " << mymath::oct_norm(A) * mymath::oct_norm(A_inv) << "\n";

		std::cout << "\n";

		std::cout << "|| A ||_inf = " << mymath::cube_norm(A) << "\n";
		std::cout << "|| A ^ -1 ||_inf = " << mymath::cube_norm(A_inv) << "\n";
		std::cout << "cond_inf A = " << mymath::cube_norm(A) * mymath::cube_norm(A_inv) << "\n";


	}
	else std::cout << "nope\n";


}

int main() {

	mymath::dynamic_matrix<test_T> A = { 
		{15, 2, -3, 7},
		{-5, 11, 2, -3},
		{0, -1, 7, 4},
		{12, 0, -6, 20}
	};

	mymath::dynamic_matrix<test_T> A_inv(0, A.rows(), A.columns());

	mymath::dynamic_vector<test_T> b = {53, -90, 107, 68};

	mymath::dynamic_matrix<test_T> B;

	mymath::dynamic_vector<test_T> x;



	try {
		x = mymath::relax_iteration(A, b, EPS);
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what() << "\n\n";
	}
	catch (std::bad_alloc& e) {
		std::cerr << "Caught bad_alloc: " << e.what() << "\n\n";
	}

	std::cout << "Matrix A: \n";
	mymath::utilities::print(A);

	std::cout << "Vector b: \n";
	mymath::utilities::print(b);
	std::cout << "\n";

	if (x.size() != 0) {

		auto y = mymath::multiply(A, x);
		std::cout << "Vector x: \n";
		mymath::utilities::print(x);
		std::cout << "\n";

		std::cout << "Vector y = Ax: \n";
		mymath::utilities::print(y);
		std::cout << "\n";

		std::cout << "|| b - y || = " << mymath::norm(b - y);
		std::cout << "\n\n\n";

		mymath::dynamic_vector<test_T> line(0, A.rows());

		line[0] = 1;
		for (size_t i = 0; i != line.size() - 1; ++i) {
			auto tmp = mymath::relax_iteration(A, line, 1e-15);

			for (size_t j = 0; j != line.size(); ++j)
				A_inv[j][i] = tmp[j];

			std::swap(line[i], line[i + 1]);
		}


		auto tmp = mymath::relax_iteration(A, line, EPS);
		for (size_t j = 0; j != line.size(); ++j)
			A_inv[j][line.size() - 1] = tmp[j];

		std::cout << "Inverse A: \n";
		mymath::utilities::print(A_inv);
		std::cout << "\n";

		std::cout << "|| A ||_1 = " << mymath::oct_norm(A) << "\n";
		std::cout << "|| A ^ -1 ||_1 = " << mymath::oct_norm(A_inv) << "\n";
		std::cout << "cond_1 A = " << mymath::oct_norm(A) * mymath::oct_norm(A_inv) << "\n";

		std::cout << "\n";

		std::cout << "|| A ||_inf = " << mymath::cube_norm(A) << "\n";
		std::cout << "|| A ^ -1 ||_inf = " << mymath::cube_norm(A_inv) << "\n";
		std::cout << "cond_inf A = " << mymath::cube_norm(A) * mymath::cube_norm(A_inv) << "\n";


	}
	else std::cout << "nope\n";



	return 0;
}
