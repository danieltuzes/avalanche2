#pragma once

#define EPSILON 1e-10 // overshoot external stress then the minimal required increment. This is not the smallest number, it is defined in gen_math.h as #define SMALLEST_NUMBER nextafter(0.,DINFTY)
typedef double deftype; // once maybe it willl be a good idea ti increase the precision to long double