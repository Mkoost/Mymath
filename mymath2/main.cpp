#include "mymath/mymath.h"
#include <iostream>

template<class T, size_t n, size_t m>
void print(const mymath::matrix<T, n, m>& A) {
	for (int i = 0; i != n; ++i) {
		for (int j = 0; j != m; ++j)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}
}

template<class T>
void print(const mymath::quaternion<T>& q) {
	std::cout << q.w << " " << q.x << " " << q.y << " " << q.z;
}


int main() {
	const mymath::imat3 B(0);
	
	mymath::imat3::fill(B);

	mymath::iquat q{1, 1, 1, 1};

	print(q);

	return 0;
}