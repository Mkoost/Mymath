#pragma once

#if defined(MYMATH_MATRIX_STATE)
#include "../headers/matrix.h"
#include "../details/__expr.inl"
#	if	MYMATH_MATRIX_STATE == 1
#		include "quaternion/math_functions_default.inl"
#		include "matrix/class_methods_default.inl"
#		include "matrix/external_operators_default.inl"
#		include "matrix/internal_classes_methods_default.inl"
#	else
#error Unsupported MYMATH_MATRIX_STATE state
#	endif
#else
#error Unsupported MYMATH_MATRIX_STATE state
#endif