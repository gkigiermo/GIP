#!/bin/sh
rm *.o
INC="-I/usr/local/cuda/include/"
INCLIB="-L/usr/local/cuda/lib64/"
mpic++ -c src/BasicInterfaces/GIP_Topo.cpp -fopenmp
mpic++ -c src/AbstractInterfaces/GIP_Vector.cpp -fopenmp
mpic++ -c src/Specifications/Cuda/GIP_VectorCuda.cpp -fopenmp $INC
mpic++ -c src/AbstractInterfaces/GIP_Matrix.cpp -fopenmp
mpic++ -c src/Specifications/Cuda/GIP_MatrixCuda.cpp -fopenmp  $INC 
mpic++ -c src/AbstractInterfaces/GIP_Timer.cpp -fopenmp
mpic++ -c src/Specifications/Cuda/GIP_TimerCuda.cpp -fopenmp  $INC
mpic++ -c src/BasicInterfaces/GIP_Parameters.cpp -fopenmp
mpic++ -c src/AbstractInterfaces/GIP_Arch.cpp -fopenmp
mpic++ -c src/Specifications/Cuda/GIP_ArchGPU.cpp -fopenmp $INC
mpic++ -c src/AbstractInterfaces/GIP_NonLinear.cpp -fopenmp
mpic++ -c src/Specifications/Cuda/GIP_NonLinearCuda.cpp -fopenmp  $INC
mpic++ -c DrivenCavity.cpp  -fopenmp $INC -DENABLE_GPU=1
mpic++ DrivenCavity.o GIP_NonLinearCuda.o GIP_NonLinear.o GIP_ArchGPU.o GIP_Arch.o GIP_Parameters.o GIP_TimerCuda.o GIP_Timer.o GIP_MatrixCuda.o GIP_Matrix.o GIP_VectorCuda.o GIP_Vector.o GIP_Topo.o -L./src/Specifications/Cuda/Libs -lgpu $INCLIB -lcudart  -o DrivenCavityGPU -fopenmp $INC
rm *.o
