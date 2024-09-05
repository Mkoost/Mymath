#pragma once

#define __MYMATH_TEMPLATE_GEN_TYPE_PSEUDONYMS_17845(root, postfix, type, ...) using d##root##postfix = type<double, __VA_ARGS__>; \
using f##root##postfix = type<float, __VA_ARGS__>; \
using i##root##postfix = type<int, __VA_ARGS__>; \
using l##root##postfix = type<long, __VA_ARGS__>; \
using ll##root##postfix = type<long long, __VA_ARGS__>; \
using u##root##postfix = type<unsigned, __VA_ARGS__>; \
using ui##root##postfix = type<unsigned int, __VA_ARGS__>; \
using ul##root##postfix = type<unsigned long, __VA_ARGS__>; \
using ull##root##postfix = type<unsigned long long, __VA_ARGS__>;

#define __GEN_MYMATH_MAT_PSEUDONYMS_10931(postfix, n, m) __MYMATH_TEMPLATE_GEN_TYPE_PSEUDONYMS_17845(mat, postfix, matrix, n, m)
#define __GEN_MYMATH_VEC_PSEUDONYMS_09182(postfix, n) __MYMATH_TEMPLATE_GEN_TYPE_PSEUDONYMS_17845(vec, postfix, vector, n)
