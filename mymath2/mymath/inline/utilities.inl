#pragma once

#include <string>
#include <fstream>

#include "matrix.inl.h"
#include "vector.inl.h"
#include "dynamic_matrix.inl.h"
#include "dynamic_vector.inl.h"
#include "../headers/data_structs.h"


namespace mymath {
	namespace utilities {
#if defined(MYMATH_MATRIX_STATE)
		template<class T, size_t n, size_t m>
		void print(const mymath::matrix<T, n, m>& A) {
			for (int i = 0; i != n; ++i) {
				for (int j = 0; j != m; ++j)
					std::cout << A[i][j] << " ";
				std::cout << "\n";
			}
			std::cout << "\n";
		}

#endif

#if defined(MYMATH_QUATERNION_STATE)
		template<class T>
		void print(const mymath::quaternion<T>& q) {
			std::cout << q.w << " " << q.x << " " << q.y << " " << q.z;
		}
#endif

#if defined(MYMATH_VECTOR_STATE)
		template<class T, size_t n>
		void print(const mymath::vector<T, n>& vec) {
			for (int i = 0; i != n; ++i)
				std::cout << vec[i] << " ";
			std::cout << "\n";
		}

		template<class T>
		void input_vector() {

		}
#endif

#if defined(MYMATH_DYNAMIC_MATRIX_STATE)
		template<class T>
		void print(const mymath::dynamic_matrix<T>& mat) {
			for (size_t i = 0, szr = mat.rows(); i != szr; ++i){
				for (size_t j = 0, szc = mat.columns(); j != szc; ++j)
					std::cout << mat[i][j] << " ";
				std::cout << "\n";
			}
			std::cout << "\n";
		}

		template<typename data_T = double>
		dynamic_matrix<data_T> input_matrix_NxN(const std::string& file) {
			std::ifstream in;

			in.open(file, std::ios::in);

			data_structs::base_data_ptr<data_T> data;
			size_t n;
			in >> n;

			data.size = n * n;

			data.ptr = new data_T[data.size];

			size_t c = 0;

			while (c < data.size && !in.eof())
				in >> data.ptr[c++];

			in.close();
			
			dynamic_matrix<data_T> lll;
			lll.move(data.ptr, n, n);
			return lll;
		}

#endif

#if defined(MYMATH_DYNAMIC_VECTOR_STATE)
		template<class T>
		void print(const mymath::dynamic_vector<T>& vec) {
			for (size_t i = 0, k = vec.size(); i != k; ++i)
				std::cout << vec[i] << " ";
			std::cout << "\n";
		}

		template<typename data_T = double>
		dynamic_vector<data_T> input_vector_N(const std::string& file) {
			std::ifstream in;

			in.open(file, std::ios::in);

			data_structs::base_data_ptr<data_T> data;
			size_t n;
			in >> n;

			data.size = n;

			data.ptr = new data_T[data.size];

			size_t c = 0;

			while (c < data.size && !in.eof())
				in >> data.ptr[c++];

			in.close();

			mymath::dynamic_vector<data_T> llll;

			llll.move(data.ptr, n);

			return llll;
		}

		
#endif


	}
}
