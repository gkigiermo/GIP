#include "GIP_ArchGPU.h"

GIP_ArchGPU::GIP_ArchGPU(int _numcpus,int _numcores,int _numdevs):GIP_Arch(_numcpus,_numcores,_numdevs)
{
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        devId=rank%_numcpus; // De momento se asume que el numero de gpus es igual la numero de cpus
}

void GIP_ArchGPU::postConstruct(int _numcpus,int _numcores, int _numdevs)
{

       int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        num_cpus=_numcpus;
        num_cores=_numcores;
        num_devs=_numdevs;
   
        devId=rank%_numcpus; // De momento se asume que el numero de gpus es igual la numero de cpus

}

void GIP_ArchGPU::setUp()
{
        GIP_cudaSetDevice(devId);
}
