#include "mymath/mymath.h"
#include <iostream>
#include <memory>
#include <chrono>     


using namespace std;

template<typename T, size_t n, size_t m>
void print(const mymath::matrix<T, n, m>& mat) {
	for (int i = 0; i != n; ++i){
		for (int j = 0; j != m; ++j)
			cout << mat.values[i][j] << " ";
		cout << endl;
	}
}

using stdmat = mymath::matrix<int, 10, 10>;
int main() {
	unique_ptr<stdmat> m{ new stdmat };
	stdmat::fill(*m, 1);


	*m += *m + (*m + *m);
	
	return 0;
}