#pragma once

#include "../../headers/quaternion.h"

namespace mymath {
	// ------------------------------------------------------------------>	 EXTERNAL OPERATIONS

	// ------------------ QUATERTION PLUS ------------------

	template<typename T>
	constexpr quaternion<T> operator+ (quaternion<T> ql, quaternion<T> qr) {
		ql += qr;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator+ (quaternion<T> ql, T val) {
		quaternion<T> tmp{ val, 0, 0, 0 }; ql += tmp;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator+ (T val, quaternion<T> qr) {
		quaternion<T> tmp{ val, 0, 0, 0 }; tmp += qr;
		return tmp;
	}

	// ------------------ QUATERTION MINUS ------------------

	template<typename T>
	constexpr quaternion<T> operator- (quaternion<T> ql, quaternion<T> qr) {
		ql -= qr;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator- (quaternion<T> ql, T val) {
		quaternion<T> tmp{ val, 0, 0, 0 }; ql -= tmp;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator- (T val, quaternion<T> qr) {
		quaternion<T> tmp{ val, 0, 0, 0 }; tmp -= qr;
		return tmp;
	}

	// ------------------ QUATERTION MUL ------------------

	template<typename T>
	constexpr quaternion<T> operator* (quaternion<T> ql, quaternion<T> qr) {
		ql *= qr;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator* (quaternion<T> ql, T val) {
		ql.w *= val;
		ql.x *= val;
		ql.y *= val;
		ql.z *= val;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator* (T val, quaternion<T> qr) {
		qr.w *= val;
		qr.x *= val;
		qr.y *= val;
		qr.z *= val;
		return qr;
	}
	// ------------------ QUATERTION DIV ------------------

	template<typename T>
	constexpr quaternion<T> operator/ (quaternion<T> ql, quaternion<T> qr) {
		ql /= qr;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator/ (quaternion<T> ql, T val) {
		ql.w /= val;
		ql.x /= val;
		ql.y /= val;
		ql.z /= val;
		return ql;
	}

	template<typename T>
	constexpr quaternion<T> operator/ (T val, quaternion<T> qr) {
		return inv(qr) * val;
	}

}