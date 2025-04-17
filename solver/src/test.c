#include <stdlib.h>
#include "fp.h"
double* getFormFactor(FP** lines, int count, int samples) {
  double* out = (double*)calloc(count, sizeof(double));
  (void)count;
  (void)samples;
  (void)lines;
  return out;
}
