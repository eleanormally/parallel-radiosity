#include <mpi.h>
#include <cassert>
#include "clockcycle.h"
#include "fp.h"
#include "mesh/argparser.h"
#include "mesh/face.h"
#include "mesh/mesh.h"
#include "mesh/meshdata.h"
#include "mesh/ray.h"
#include "mesh/raytracer.h"
#include "mesh/vectors.h"

#define sample_points 1024
#define shadow_samples 32
#define EPISLON 0.00001
using namespace std;

f3 vecToF3(const Vec3f& v) {
  f3 out;
  out.data[0] = v.x();
  out.data[1] = v.y();
  out.data[2] = v.z();
  return out;
}
double* computeFormFactors(int i, ArgParser& args, int rank) {
  Mesh* m = args.mesh;
  FP**  lines = (FP**)calloc(m->numFaces(), sizeof(FP*));

  Face* f = m->getFace(i);

  for (int j = 0; j < m->numFaces(); j++) {
    lines[j] = (FP*)calloc(sample_points, sizeof(FP));

    Face* f2 = m->getFace(j);
    f3    anorm = vecToF3(f->computeNormal());
    f3    bnorm = vecToF3(f2->computeNormal());
    for (int k = 0; k < sample_points; k++) {
      f3 a = vecToF3(f->RandomPoint());
      f3 b = vecToF3(f2->RandomPoint());
      lines[j][k] = FP(a, anorm, b, bnorm);
    }
  }
  double* output = getFormFactor(lines, m->numFaces(), sample_points, rank);
  for (int j = 0; j < m->numFaces(); j++) {
    output[j] /= f->getArea();
    free(lines[j]);
  }
  for (int j = 0; j < m->numFaces(); j++) {
    Face*  f2 = m->getFace(j);
    double totalVisible = 0.0;
    for (int k = 0; k < shadow_samples; k++) {
      Vec3f iPoint = f->RandomPoint();
      Vec3f jPoint = f2->RandomPoint();
      Vec3f iToJ = jPoint - iPoint;

      const double iToJLength = iToJ.Length();
      iToJ.Normalize();
      const Ray r(iPoint, iToJ);
      Hit       h;
      args.raytracer->CastRay(r, h, true);
      if (fabs(h.getT() - iToJLength) < EPSILON) {
        totalVisible += 1.0;
      }
    }
    totalVisible /= shadow_samples;
    output[j] *= totalVisible;
  }

  double flux = 0;
  for (int j = 0; j < m->numFaces(); j++) {
    flux += output[j];
  }
  for (int j = 0; j < m->numFaces(); j++) {
    output[j] /= flux;
  }
  free(lines);
  return output;
}

int main(int argc, const char* argv[]) {
  MPI_File fh;
  MPI_Init(nullptr, nullptr);
  int worldRank;
  int worldSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Barrier(MPI_COMM_WORLD);
  MeshData  mesh{};
  ArgParser args(argc, argv, &mesh);
  Mesh*     m = args.mesh;
  for (int i = 0; i < args.num_subdivisions; i++) {
    m->Subdivision();
  }

  uint64_t start = clock_now();

  //creating output file
  MPI_File_open(MPI_COMM_WORLD, "formfactors.out",
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  if (worldRank == 0) {
    int faceCount = m->numFaces();
    MPI_File_write_at(fh, 0, &faceCount, 1, MPI_INT, MPI_STATUS_IGNORE);
  }
  uint64_t totalWriteTime = 0;

  for (int i = 0; i < m->numFaces(); i++) {
    if (i % worldSize == worldRank) {
      double* out = computeFormFactors(i, args, worldRank);
      for (int j = 0; j < m->numFaces(); j++) {
        int offset = sizeof(int) + ((i * m->numFaces() + j) * sizeof(double));
        uint64_t write = clock_now();
        MPI_File_write_at(fh, offset, &out[j], sizeof(double), MPI_DOUBLE,
                          MPI_STATUS_IGNORE);
        uint64_t endwrite = clock_now();
        totalWriteTime += (endwrite - write);
      }
      free(out);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_close(&fh);
  uint64_t endTime = clock_now();
  if (worldRank == 0) {
    printf("total time: %llu, non-IO Time: %llu, IO Time: %llu\n",
           endTime - start, endTime - start - totalWriteTime, totalWriteTime);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
