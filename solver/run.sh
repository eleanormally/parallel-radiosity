#!/bin/bash -x

#SBATCH --partition=el8-rpi
#SBATCH --gres=gpu:4
#SBATCH --nodes=4
#SBATCH --time=120

cd ~/scratch/build/
module load spectrum-mpi cuda/11.2
make
echo "Weak Scaling"
mpirun -n 1 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 0
mpirun -n 4 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 1
mpirun -n 16 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 2
echo "Strong Scaling"
mpirun -n 1 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 3
mpirun -n 2 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 3
mpirun -n 4 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 3
mpirun -n 8 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 3
mpirun -n 16 ./radiosity --input ../inputs/cornell_box_diffuse_sphere.obj --subdiv 3
