#include<stdio.h>
#include<iostream>
#include "src/Algorithms/GIP_DNS.h"

#if ENABLE_GPU 
    #include "src/Specifications/Cuda/Cuda.h"
    typedef GIP_MatrixCuda Matrix;
    typedef GIP_VectorCuda Vector;
    typedef GIP_NonLinearCuda NonLinear;
    typedef GIP_ArchGPU Architecture;
#else
    #include "src/Specifications/Multicore/Multicore.h"
    typedef GIP_MatrixMulticore Matrix;
    typedef GIP_VectorMulticore Vector;
    typedef GIP_NonLinearMulticore NonLinear;
    typedef GIP_ArchCPU Architecture;

#endif

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

    GIP_DNS<Matrix,Vector,Architecture,NonLinear> driven(name);

    driven.setUp();
 
    driven.integrate();   
    

    MPI_Finalize();
}
