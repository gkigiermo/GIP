#!/bin/sh
mpic++ -c src/BasicInterfaces/GIP_Topo.cpp -fopenmp
mpic++ -c src/AbstractInterfaces/GIP_Vector.cpp -fopenmp
mpic++ -c src/Specifications/Multicore/GIP_VectorMulticore.cpp -fopenmp
mpic++ -c src/AbstractInterfaces/GIP_Matrix.cpp -fopenmp
mpic++ -c src/Specifications/Multicore/GIP_MatrixMulticore.cpp -fopenmp
mpic++ -c src/AbstractInterfaces/GIP_Timer.cpp -fopenmp
mpic++ -c src/Specifications/Multicore/GIP_TimerMulticore.cpp -fopenmp
mpic++ -c src/BasicInterfaces/GIP_Parameters.cpp -fopenmp
mpic++ -c src/AbstractInterfaces/GIP_Arch.cpp -fopenmp
mpic++ -c src/Specifications/Multicore/GIP_ArchCPU.cpp -fopenmp
mpic++ -c DrivenCavity.cpp  -fopenmp
mpic++ DrivenCavity.o GIP_ArchCPU.o GIP_Arch.o GIP_Parameters.o GIP_TimerMulticore.o GIP_Timer.o GIP_MatrixMulticore.o GIP_Matrix.o GIP_VectorMulticore.o GIP_Vector.o GIP_Topo.o -o DrivenCavityCPU.x -fopenmp
rm *.o
