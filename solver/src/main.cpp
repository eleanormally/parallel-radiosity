#include <mpi.h>
#include <cassert>
#include <fstream>
#include "mesh/vectors.h"

int main(int argc, char* argv[]) {
  assert(argc >= 2);
  std::ifstream file;
  file.open(argv[1]);
  std::string line;

  int dimension;
  file >> dimension;
  std::vector<std::vector<double>> formFactors(
      dimension, std::vector<double>(dimension, 0.0));
  std::vector<Vec3f> emmitence(dimension);
  std::vector<Vec3f> reflectance(dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      file >> formFactors[i][j];
    }
    int x, y, z, r, g, b;
    file >> x >> y >> z >> r >> g >> b;
    emmitence[i] = Vec3f(x, y, z);
    reflectance[i] = Vec3f(r, g, b);
  }
  MPI_Init(nullptr, nullptr);
  int worldRank;
  int worldSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Barrier(MPI_COMM_WORLD);
}
