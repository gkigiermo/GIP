#!/bin/sh
rm *.o
mpic++ -c BasicInterfaces/GIP_Topo.cpp -fopenmp
mpic++ -c AbstractInterfaces/GIP_Vector.cpp -fopenmp
mpic++ -c Specifications/Cuda/GIP_VectorCuda.cpp -fopenmp
mpic++ -c AbstractInterfaces/GIP_Matrix.cpp -fopenmp
mpic++ -c Specifications/Cuda/GIP_MatrixCuda.cpp -fopenmp
mpic++ -c AbstractInterfaces/GIP_Timer.cpp -fopenmp
mpic++ -c Specifications/Cuda/GIP_TimerCuda.cpp -fopenmp
mpic++ -c BasicInterfaces/GIP_Parameters.cpp -fopenmp
mpic++ -c AbstractInterfaces/GIP_Arch.cpp -fopenmp
mpic++ -c Specifications/Cuda/GIP_ArchGPU.cpp -fopenmp
mpic++ -c DrivenCavity.cpp  -fopenmp
mpic++ DrivenCavity.o GIP_ArchGPU.o GIP_Arch.o GIP_Parameters.o GIP_TimerCuda.o GIP_Timer.o GIP_MatrixCuda.o GIP_Matrix.o GIP_VectorCuda.o GIP_Vector.o GIP_Topo.o -L./Specifications/Cuda/Libs -lgpu -lcudart  -o DrivenCavityGPU.x -fopenmp
rm *.o
