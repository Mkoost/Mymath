#include "mymath/mymath.h"
#include <iostream>


using test_T = double;
constexpr const size_t SIZE = 4;
constexpr const double EPS = 1e-12;


int main(){
	// Eigenvalues var 9 {107.973, -38.7335, 32.1177, 14.4424}
	//mymath::utilities::input_matrix_NxN("Test/lab 3/matrix.dat");

	mymath::dynamic_matrix<test_T> A = { 
		{1.5, 0, -0.43, -0.75},
		{0, 3, 0.87, -0.5}, 
		{-0.43, 0.87, 2.9, -0.22}, 
		{-0.75, -0.5, -0.22, 2.6}}; 

	mymath::utilities::print(A);
	std::cout << "\n";
	mymath::francis_eigenvals(A);

	mymath::utilities::print(A);

	return 0;
}