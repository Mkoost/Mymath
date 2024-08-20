#include "mymath/mymath.h"
#include <iostream>
#include <memory>
#include <chrono>     


using namespace std;


using quat = mymath::quaternion<int>;
using mat3 = mymath::matrix<int, 3, 3>;

int main() {
	mat3 m;
	mat3::eve(m);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
			cout << m.values[i][j] << " ";
		cout << endl;
	}
	return 0;
}