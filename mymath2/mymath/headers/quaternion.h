#pragma once

namespace mymath {

	template<typename T>
	struct quaternion {
		T w, x, y, z;

		constexpr quaternion<T>& operator+= (const quaternion<T>& q);

		constexpr quaternion<T>& operator-= (const quaternion<T>& q);

		constexpr quaternion<T>& operator*= (const quaternion<T>& q);

		constexpr quaternion<T>& operator/= (const quaternion<T>& q);
	};

	// -------------------------------------------------------->		 MATH FUNCTIONS / OPERATIONS  

	constexpr const quaternion<double> I{ 0, 1, 0, 0 };
	constexpr const quaternion<double> J{ 0, 0, 1, 0 };
	constexpr const quaternion<double> K{ 0, 0, 0, 1 }; //noexcept

	using dquat = quaternion<double>;
	using fquat = quaternion<float>;

	using iquat = quaternion<int>; 
	using lquat = quaternion<long>;
	using llquat = quaternion<long long>; 

	using uquat = quaternion<unsigned>;
	using uiquat = quaternion<unsigned int>;
	using ulquat = quaternion<unsigned long>;
	using ullquat = quaternion<unsigned long long>;

	using quat = quaternion<int>;
}

