#pragma once

#define MYMATH_INCLUDE_SETTINGS

// Header includes if state equals 1
#define MYMATH_QUATERNION_STATE 1
#define MYMATH_MATRIX_STATE 1
#define MYMATH_VECTOR_STATE 1

// The matrix is multiplied using:
// 1) the Strassen algorithm if state equals 1
// 2) Standart algorithm if state equals 2

#define MYMATH_MATRIX_S_MUL 1 


