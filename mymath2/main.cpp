#include "mymath/mymath.h"
#include <iostream>


using test_T = double;
constexpr const size_t SIZE = 4;

void qr_solve() {
	mymath::matrix<test_T, SIZE, SIZE> A{	
		{ 28.859, -0.008, 2.406, 19.240 },
		{ 14.436, -0.001, 1.203, 9.624 },
		{ 120.204, -0.032, 10.024, 80.144 },
		{ -57.714, 0.016, -4.812, -38.478 }};

	mymath::matrix<test_T, SIZE, SIZE> B;

	mymath::vector<test_T, SIZE> b = { 30.459, 18.248, 10, -60.908 };
	mymath::vector<test_T, SIZE>* x = nullptr;

	try {
		x = mymath::qr_solve <test_T, SIZE, test_T>(A, b, nullptr, &B, 1e-9);
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

	if (x != nullptr) {

		auto y = mymath::multiply(A, *x);
		std::cout << "Vector x: \n";
		mymath::utilities::print(*x);
		std::cout << "\n";

		std::cout << "Vector Ax: \n";
		mymath::utilities::print(*y);
		std::cout << "\n";

		std::cout << "Norm b - y: \n";
		std::cout << mymath::norm(b - *y);
		std::cout << "\n";

		delete x;
		delete y;
	}
	else std::cout << "nullptr\n";
}

int main() {
	qr_solve();
	return 0;
}