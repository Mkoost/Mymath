#pragma once

#if defined(MYMATH_DYNAMIC_VECTOR_STATE)
#include "../headers/dynamic_vector.h"
#include "../details/__expr.inl"
#	if	MYMATH_DYNAMIC_VECTOR_STATE == 1
#		include "dynamic_vector/internal_classes_methods_default.inl"
#		include "dynamic_vector/math_functions_default.inl"
#		include "dynamic_vector/class_methods_default.inl"
#		include "dynamic_vector/external_operators_default.inl"
#	else
#error Unsupported MYMATH_DYNAMIC_VECTOR_STATE state
#	endif
#endif