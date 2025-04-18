#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include "fp.h"

f3 inline __device__ sub(f3 a, f3 b) {
  f3 out;
  out.data[0] = a.data[0]-b.data[0];
  out.data[1] = a.data[1]-b.data[1];
  out.data[2] = a.data[2]-b.data[2];
  return out;
}
double inline __device__ dot(f3 a, f3 b) {
  double out = 0;
  out += a.data[0]*b.data[0];
  out += a.data[1]*b.data[1];
  out += a.data[2]*b.data[2];
  return out;
}
double inline __device__ len(f3 a) {
  double len = 0;
  len += a.data[0]*a.data[0];
  len += a.data[1]*a.data[1];
  len += a.data[2]*a.data[2];
  return sqrt(len);
}



void __global__ calculate(FP** lines, double* out) {
  extern __shared__ double lineFactor[];
  FP p = lines[blockIdx.x][threadIdx.x];
  f3 iToJ = sub(p.a,p.b);
  f3 jToI = sub(p.b,p.a);
  double iCos = dot(p.anorm, iToJ) / len(iToJ);
  double jCos = dot(p.bnorm, jToI) / len(jToI);
  lineFactor[threadIdx.x] = fabs(iCos*jCos)/ (len(iToJ)*len(jToI));
  __syncthreads();
  for(unsigned int s=blockDim.x/2; s > 0; s >>= 1) {
    if(threadIdx.x < s) {
      lineFactor[threadIdx.x] += lineFactor[threadIdx.x+s];
    }
    __syncthreads();
  }
  if(threadIdx.x == 0) {
    out[blockIdx.x] = lineFactor[0];
  }
}

double* getFormFactor(FP** lines, int count, int samples, int rank) {
  int c;
  cudaGetDeviceCount(&c);
  cudaSetDevice(rank % c);
  double* cudaOut;
  FP** cudaLines;
  cudaMallocManaged(&cudaLines, count*sizeof(FP*));
  for(int i = 0; i < count; i++) {
    cudaMallocManaged(cudaLines+i,samples*sizeof(FP));
    for(int j = 0; j < samples; j++) {
      cudaLines[i][j] = lines[i][j];
    }
  }
  cudaMallocManaged(&cudaOut, count*sizeof(double));
  calculate<<<count, samples, samples*sizeof(double)>>>(cudaLines, cudaOut);
  cudaDeviceSynchronize();
  double* out = (double*)calloc(count, sizeof(double));
  for(int i = 0; i < count; i++) {
    out[i] = cudaOut[i];
    cudaFree(cudaLines[i]);
  }
  cudaFree(cudaLines);
  cudaFree(cudaOut);

  (void)count;
  (void)samples;
  (void)lines;
  return out;
}
