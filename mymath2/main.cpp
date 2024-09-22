#include "mymath/mymath.h"
#include <iostream>


using test_T = double;
constexpr const size_t SIZE = 4;


template<class T, class tmp_T>
using solver_ptr = mymath::data_structs::base_data_dynamic_vector_matrix<T, 1, 1>(*)(mymath::dynamic_matrix<T>&, mymath::dynamic_vector<T>&, tmp_T);

template<class T, class tmp_T>
mymath::data_structs::base_data_dynamic_vector_matrix<T, 1, 1> s(mymath::dynamic_matrix<T>& A, mymath::dynamic_vector<T>& b, tmp_T zero) {
	return mymath::gauss_solve(A, b, zero);
};

template<class T, class tmp_T>
void test_solve(const std::string& mat_path, const std::string& vec_path, solver_ptr<T, tmp_T> solver) {
	mymath::dynamic_matrix<T> A = mymath::utilities::input_matrix_NxN<T>(mat_path);

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

		std::cout << "|| b - y ||: \n";
		std::cout << mymath::norm(b - y);
		std::cout << "\n";

	}
	else std::cout << "nope\n";


}

int main() {
	std::string path1 = "C:\\Users\\Миша\\source\\repos\\mymath2\\mymath2\\Test\\lab1\\sys5\\matrix.dat";
	std::string path2 = "C:\\Users\\Миша\\source\\repos\\mymath2\\mymath2\\Test\\lab1\\sys5\\vector.dat";


	test_solve<test_T, test_T>(path1, path2, s);

	return 0;
}