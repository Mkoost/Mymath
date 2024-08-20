#ifndef MYMATH_QUATERNION
#define MYMATH_QUATERNION
#include <cmath>

namespace mymath {

	template<typename T>
	struct quaternion {
		T w, x, y, z;

		quaternion<T>& operator+= (const quaternion<T>& q);

		quaternion<T>& operator-= (const quaternion<T>& q);

		quaternion<T>& operator*= (const quaternion<T>& q);

		quaternion<T>& operator/= (const quaternion<T>& q);
	};

	// -------------------------------------------------------->		 MATH FUNCTIONS / OPERATIONS  

	const quaternion<double> I{ 0, 1, 0, 0 };
	const quaternion<double> J{ 0, 0, 1, 0 };
	const quaternion<double> K{ 0, 0, 0, 1 };

}

#endif
