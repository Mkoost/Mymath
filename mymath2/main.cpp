#include "mymath/mymath.h"
#include <iostream>


int main() {	
	mymath::imat3 A{
		{1, 6, 1},
		{1, 1, 3},
		{3, 1, 3} };
	mymath::ivec3 v{1, 1, 1};

	mymath::ivec3 c;

	mymath::multiply(v, A, &c);

	mymath::utilities::print(c);


	return 0;
}