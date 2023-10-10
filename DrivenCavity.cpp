#include<stdio.h>
#include<iostream>
#include "src/Specifications/Multicore/Multicore.h"
//#include "src/Specifications/Cuda/Cuda.h"


#include "src/Algorithms/GIP_DNS.h"


//typedef GIP_MatrixCuda Matrix;
//typedef GIP_VectorCuda Vector;
//typedef GIP_ArchGPU Architecture;

typedef GIP_MatrixMulticore Matrix;
typedef GIP_VectorMulticore Vector;
typedef GIP_ArchCPU Architecture;


using namespace std;
int main(int argc, char** argv)
{

    if(argc-1 != 1)
    {
        cout<<"a.out <Mesh_file>  "<<endl;
        return 0;
    }

    int prov;
    int req=MPI_THREAD_SINGLE;
    MPI_Init_thread(0,0,req, &prov);
    char* name=argv[1];
    GIP_DNS<Matrix,Vector,Architecture> driven(name);
    driven.setUp();
    driven.integrate();   


    MPI_Finalize();
}
