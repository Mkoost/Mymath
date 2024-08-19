#if defined(MYMATH_QUATERNION_STATE) && (MYMATH_QUATERNION_STATE == 1)
#ifndef MYMATH_QUATERNION
#define MYMATH_QUATERNION

#include <cmath>

namespace mymath {

	// ------------------------------------ CLASSES ------------------------------------

	template<typename T>
	class quaternion {
	public:
		T w, x, y, z;

		quaternion<T>& operator+= (const quaternion<T>& q) {
			x += q.x; y += q.y; w += q.w; z += q.z;
			return *this;
		}

		quaternion<T>& operator-= (const quaternion<T>& q) {
			x -= q.x; y -= q.y; w -= q.w; z -= q.z;
			return *this;
		}
		quaternion<T>& operator*= (const quaternion<T>& q) {
			w = w * q.w + x * q.x + y * q.y + z * q.z;
			x = w * q.x + x * q.w + y * q.z - z * q.y;
			y = w * q.y + y * q.w + z * q.x - x * q.z;
			z = w * q.z + z * q.w + x * q.y - y * q.x;
			return *this;
		};

		quaternion<T>& operator/= (const quaternion<T>& q) {
			*this *= conj(q);
			T tmp = norm2(q);
			w /= tmp;
			x /= tmp;
			y /= tmp;
			z /= tmp;
			return *this;
		}

	};

	// ------------------------------------ QUATERTION CONSTS ------------------------------------

	const quaternion<double> I{ 0, 1, 0, 0 };
	const quaternion<double> J{ 0, 0, 1, 0 };
	const quaternion<double> K{ 0, 0, 0, 1 };



	// ------------------------------------ OPERATORS ------------------------------------
	

	// ------------------ QUATERTION PLUS ------------------

	template<typename T>
	quaternion<T> operator+ (quaternion<T> ql, quaternion<T> qr) {
		ql += qr;
		return ql;
	}

	template<typename T>
	quaternion<T> operator+ (quaternion<T> ql, T val) {
		quaternion<T> tmp{ val, 0, 0, 0 }; ql += tmp;
		return ql;
	}

	template<typename T>
	quaternion<T> operator+ (T val, quaternion<T> qr) {
		quaternion<T> tmp{ val, 0, 0, 0 }; tmp += qr;
		return tmp;
	}

	// ------------------ QUATERTION MINUS ------------------

	template<typename T>
	quaternion<T> operator- (quaternion<T> ql, quaternion<T> qr) {
		ql -= qr;
		return ql;
	}

	template<typename T>
	quaternion<T> operator- (quaternion<T> ql, T val) {
		quaternion<T> tmp{ val, 0, 0, 0 }; ql -= tmp;
		return ql;
	}

	template<typename T>
	quaternion<T> operator- (T val, quaternion<T> qr) {
		quaternion<T> tmp{ val, 0, 0, 0 }; tmp -= qr;
		return tmp;
	}

	// ------------------ QUATERTION MUL ------------------

	template<typename T>
	quaternion<T> operator* (quaternion<T> ql, quaternion<T> qr) {
		ql *= qr;
		return ql;
	}

	template<typename T>
	quaternion<T> operator* (quaternion<T> ql, T val) {
		ql.w *= val;
		ql.x *= val;
		ql.y *= val;
		ql.z *= val;
		return ql;
	}

	template<typename T>
	quaternion<T> operator* (T val, quaternion<T> qr) {
		qr.w *= val;
		qr.x *= val;
		qr.y *= val;
		qr.z *= val;
		return qr;
	}
	// ------------------ QUATERTION DIV ------------------

	template<typename T>
	quaternion<T> operator/ (quaternion<T> ql, quaternion<T> qr) {
		ql /= qr;
		return ql;
	}

	template<typename T>
	quaternion<T> operator/ (quaternion<T> ql, T val) {
		ql.w /= val;
		ql.x /= val;
		ql.y /= val;
		ql.z /= val;
		return ql;
	}

	template<typename T>
	quaternion<T> operator/ (T val, quaternion<T> qr) {
		return inv(qr) * val;
	}

	// ------------------------------------ MATH FUNCTIONS / OPERATIONS  ------------------------------------

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

#endif
#endif