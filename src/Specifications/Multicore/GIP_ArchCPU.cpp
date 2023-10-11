#include "GIP_ArchCPU.h"

GIP_ArchCPU::GIP_ArchCPU(int _numcpus,int _numcores,int _numdevs):GIP_Arch(_numcpus,_numcores,_numdevs)
{
    devId=0;
}

void GIP_ArchCPU::postConstruct(int _numcpus,int _numcores, int _numdevs)
{

       int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        num_cpus=_numcpus;
        num_cores=_numcores;
        num_devs=_numdevs;
   

}

void GIP_ArchCPU::setUp()
{
    omp_set_num_threads(num_cores);

    //int cores=num_cores*num_cpus; //12;
    int nthreads=num_cores;//6;
    omp_set_num_threads(nthreads);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef AFFINITY
    int cpu_id=rank%num_cpus;
#pragma omp parallel default(shared) num_threads(nthreads)
    {  

        cpu_set_t mask_SYSTEM;
        unsigned char *ptr = (unsigned char *)(&mask_SYSTEM);
        int num_th = omp_get_num_threads();
        int id_th = omp_get_thread_num();
        for(int j=0; j<num_th; j++){
            if(id_th==j){
                if(sched_getaffinity(0, sizeof(cpu_set_t), &mask_SYSTEM))
                {
                    cout<<"TrySetAffinity: sched_getaffinity fail!"<<endl;
                }
                for(int i=0; i<cores; i++){
                    if(i != id_th+ (cpu_id)*nthreads){SetBitListOffMacro(ptr, i);/* printf("0 ");*/}
                    else{ SetBitListOnMacro(ptr, i);/* printf("1 ");*/}    
                }
                if(sched_setaffinity(0, sizeof(cpu_set_t), &mask_SYSTEM))
                {
                    cout<<"TrySetAffinity: sched_setaffinity fail!"<<endl;
                }    
            }
#pragma omp barrier 
        }
    }
#endif
}
