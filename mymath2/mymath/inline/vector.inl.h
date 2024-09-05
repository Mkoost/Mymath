#pragma once


#if defined(MYMATH_VECTOR_STATE)
#	include "../headers/vector.h"
#	if	MYMATH_VECTOR_STATE == 1
#		include "vector/math_functions_default.inl"
#		include "vector/class_methods_default.inl"
#		include "vector/external_operators_default.inl"
#		include "vector/internal_classes_methods_default.inl"
#	else
#error Unsupported MYMATH_VECTOR_STATE state
#	endif
#endif