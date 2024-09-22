#pragma once

#if defined(MYMATH_DYNAMIC_MATRIX_STATE)
#include "../headers/dynamic_matrix.h"
#include "../details/__expr.inl"
#	if	MYMATH_DYNAMIC_MATRIX_STATE == 1
#		include "dynamic_matrix/internal_classes_methods_default.inl"
#		include "dynamic_matrix/math_functions_default.inl"
#		include "dynamic_matrix/class_methods_default.inl"
#		include "dynamic_matrix/external_operators_default.inl"
#	else
#error Unsupported MYMATH_DYNAMIC_MATRIX_STATE state
#	endif
#endif