#include <mpi.h>
#include <cassert>
#include "mesh/argparser.h"
#include "mesh/mesh.h"

int main(int argc, char* argv[]) {
  assert(argc >= 2);
  Mesh mesh;
  mesh.Load(argv[1]);
  return 0;
  MPI_Init(nullptr, nullptr);
  int worldRank;
  int worldSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Barrier(MPI_COMM_WORLD);
}
