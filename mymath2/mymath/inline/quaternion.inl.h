#pragma once


#if defined(MYMATH_QUATERNION_STATE)
#	include "../headers/quaternion.h"
#	if	MYMATH_QUATERNION_STATE == 1
#		include "quaternion/math_functions_default.inl"
#		include "quaternion/class_methods_default.inl"
#		include "quaternion/external_operators_default.inl"
#	else
#error Unsupported MYMATH_QUATERNION_STATE state
#	endif
#endif