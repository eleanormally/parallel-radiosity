#ifndef __fp_h
#define __fp_h
#include <cmath>
typedef struct f3 {
  double data[3];
} f3;
typedef struct FP {
  f3 a;
  f3 anorm;
  f3 b;
  f3 bnorm;
  FP(f3 _a, f3 _anorm, f3 _b, f3 _bnorm)
      : a(_a), anorm(_anorm), b(_b), bnorm(_bnorm) {}

} FP;
double* getFormFactor(FP** lines, int count, int samples, int rank);
#endif
