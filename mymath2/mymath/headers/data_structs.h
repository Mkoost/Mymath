#pragma once

#include "./dynamic_matrix.h"
#include "./dynamic_vector.h"

namespace mymath {
	namespace data_structs{
		template<typename data_T>
		struct base_data_ptr {
			data_T* ptr = nullptr;
			size_t size = 0;
		};

		template<typename data_T, size_t N = 1, size_t K = 1>
		struct base_data_dynamic_vector_matrix {
			dynamic_matrix<data_T> mat[N];
			dynamic_vector<data_T> vec[K];
		};
	
	}
}
