#pragma once

#include "../../headers/quaternion.h"

namespace mymath{

	// ------------------------------------------------------------------>	 CLASS METHODS 

	template<typename T>
	constexpr quaternion<T>& quaternion<T>::operator+= (const quaternion<T>& q) {
		x += q.x; y += q.y; w += q.w; z += q.z;
		return *this;
	};

	template<typename T>
	constexpr quaternion<T>& quaternion<T>::operator-= (const quaternion<T>& q) {
		x -= q.x; y -= q.y; w -= q.w; z -= q.z;
		return *this;
	};

	template<typename T>
	constexpr quaternion<T>& quaternion<T>::operator*= (const quaternion<T>& q) {
		w = w * q.w + x * q.x + y * q.y + z * q.z;
		x = w * q.x + x * q.w + y * q.z - z * q.y;
		y = w * q.y + y * q.w + z * q.x - x * q.z;
		z = w * q.z + z * q.w + x * q.y - y * q.x;
		return *this;
	};

	template<typename T>
	constexpr quaternion<T>& quaternion<T>::operator/= (const quaternion<T>& q) {
		*this *= conj(q);
		T tmp = norm2(q);
		w /= tmp;
		x /= tmp;
		y /= tmp;
		z /= tmp;
		return *this;
	}
}