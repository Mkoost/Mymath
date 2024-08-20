#ifndef MYMATH
#define MYMATH

#include "settings.h"


#if defined(MYMATH_QUATERNION_STATE)
#	include "headers/quaternion.h"
#	if MYMATH_QUATERNION_STATE == 1
#		include "inline/quaternion_default.inl"
#	else
#		include "inline/quaternion_default.inl"
#	endif
#endif

#if defined(MYMATH_MATRIX_STATE) && (MYMATH_MATRIX_STATE == 1)
#	include "headers/matrix.h"
#	if MYMATH_MATRIX_STATE == 1
#		include "inline/matrix_default.inl"
#	else
#		include "inline/matrix_default.inl"
#	endif
#endif



#endif