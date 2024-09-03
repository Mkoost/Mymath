#pragma once
#include "./../../headers/matrix.h"

namespace mymath {

	template<class T, size_t n, size_t m>
	T& matrix<T, n, m>::row::operator[](size_t i) {
		return ptr[i];
	}

	template<class T, size_t n, size_t m>
	T matrix<T, n, m>::const_row::operator[](size_t i) const {
		return ptr[i];
	}


}