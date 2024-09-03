#pragma once

#include "../../headers/quaternion.h"

namespace mymath {
	// ------------------------------------------------------------------>	 MATH FUNCTIONS / OPERATIONS

	template<typename T>
	T norm2(quaternion<T> q) {
		return q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
	}

	template<typename T>
	T norm(quaternion<T> q) {
		return std::sqrt(norm2(q));
	}

	template<typename T>
	quaternion<T> conj(quaternion<T> q) {
		return quaternion<T>{q.w, -q.x, -q.y, -q.z};
	}

	template<typename T>
	quaternion<T> inv(quaternion<T> q) {
		return conj(q) / norm2(q);
	}
}