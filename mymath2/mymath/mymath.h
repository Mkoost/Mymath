#pragma once

#include <stdexcept>
#include <cmath>
#include <memory>
#include <iostream>
#include <string>
#include <fstream>

#include "settings.h"

#include "./details/__expr.inl"
#include "./details/__macroses.inl"
#include "./headers/data_structs.h"

#include "inline/quaternion.inl.h"
#include "inline/matrix.inl.h"
#include "inline/vector.inl.h"
#include "inline/utilities.inl"
#include "inline/math_functions.inl"
#include "inline/dynamic_matrix.inl.h"
#include "inline/dynamic_vector.inl.h"

#undef __MYMATH_TEMPLATE_GEN_TYPE_PSEUDONYMS_17845
#undef __GEN_MYMATH_MAT_PSEUDONYMS_10931
#undef __GEN_MYMATH_VEC_PSEUDONYMS_09182