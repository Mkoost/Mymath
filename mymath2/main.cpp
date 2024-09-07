#include "mymath/mymath.h"
#include <iostream>


int main() {	

	using dvec1000 = mymath::vector<double, 1000>;
	auto a = new dvec1000;
	dvec1000::fill(*a);

	auto b = new dvec1000;
	dvec1000::fill(*b, 1);

	*a += *b;

	mymath::utilities::print(*a);

	delete a;
	delete b;

	return 0;
}