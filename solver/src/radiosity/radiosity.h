#include <cstddef>
#include <vector>
#include "../utils/vectors.h"
using std::size_t;
using std::vector;

typedef struct RadiosityInput {
  vector<vector<float>> formFactors;
  vector<Vec3f>         emmittence;
  vector<Vec3f>         diffuse;
} RadiosityInput;
