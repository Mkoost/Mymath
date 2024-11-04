#include "mymath/mymath.h"
#include <iostream>


using test_T = double;
constexpr const size_t SIZE = 4;
constexpr const double EPS = 1e-5;


int main(){
	// First exemple {4.01142, 2.98702, 2.00425, 0.997313}
	// Eigenvalues var 9 {107.97337296690272, -38.73345211444923, 32.11767309001762, 14.44240605752893}
	//mymath::utilities::input_matrix_NxN("Test/lab 3/matrix.dat");
	/*
	{ 
		{1.5, 0, -0.43, -0.75},
		{0, 3, 0.87, -0.5}, 
		{-0.43, 0.87, 2.9, -0.22}, 
		{-0.75, -0.5, -0.22, 2.6}};
	*/

	mymath::dynamic_matrix<test_T> A = {
		{1.5, 0, -0.43, -0.75},
		{0, 3, 0.87, -0.5},
		{-0.43, 0.87, 2.9, -0.22},
		{-0.75, -0.5, -0.22, 2.6} };


	mymath::utilities::print(A);
	std::cout << "\n";

	mymath::dynamic_vector<test_T> res = mymath::francis_eigenvals(A, EPS);

	mymath::utilities::print(res);

	mymath::dynamic_vector<test_T> ans = { 107.97337296690272, -38.73345211444923, 32.11767309001762, 14.44240605752893 };

	
	std::cout << mymath::cube_norm(ans - res) << "\n\n\n";
	
	mymath::utilities::print(A);

	mymath::dynamic_vector<test_T> x0 = { -0.87, 0.00, -0.25, -0.43 };

	x0 = mymath::inverse_iteration(A, x0, 1e-3);
	
	mymath::utilities::print(x0);


	
	

	return 0;
}