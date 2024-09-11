#include "mymath/mymath.h"
#include <iostream>


int main() {

	mymath::matrix<float, 5, 5> A{
		{1, 2, 3, 4, 5},
		{5, 1, 2, 3, 4},
		{4, 5, 1, 2, 3},
		{3, 4, 5, 1, 2},
		{2, 3, 4, 5, 1} };
	
	mymath::matrix<float, 5, 5> B{
		{1, 2, 3, 4, 5},
		{5, 1, 2, 3, 4},
		{4, 5, 1, 2, 3},
		{3, 4, 5, 1, 2},
		{2, 3, 4, 5, 1} };

	mymath::vector<float, 5> b = { 1, 1, 1, 1, 1 };
	
	auto x = mymath::qr_solve(A, b);

	mymath::utilities::print(A);
	std::cout << "\n";
	mymath::utilities::print(b);
	std::cout << "\n";
	
	
	auto y = mymath::multiply(B, *x);

	mymath::utilities::print(*x);
	std::cout << "\n";
	mymath::utilities::print(*y);
	std::cout << "\n" << mymath::norm(b - *y);

	if (x) delete x;
	if (y) delete y;
	return 0;
}