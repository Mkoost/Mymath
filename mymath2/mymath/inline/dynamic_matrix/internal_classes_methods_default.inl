#pragma once
#include "./../../headers/dynamic_matrix.h"

namespace mymath {

	template<class T>
	dynamic_matrix<T>::row::row(T* ptr_) { ptr = ptr_; }

	template<class T>
	T& dynamic_matrix<T>::row::operator[](size_t i) {
		return ptr[i];
	}

	
	template<class T>
	dynamic_matrix<T>::const_row::const_row(T const * ptr_) { ptr = ptr_; }

	template<class T>
	T dynamic_matrix<T>::const_row::operator[](size_t i) const {
		return ptr[i];
	}


}