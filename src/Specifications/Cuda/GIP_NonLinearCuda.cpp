#include "GIP_NonLinearCuda.h"

double GIP_NonLinearCuda::getCFL(GIP_Topo* myTopo,GIP_Vector* u,GIP_Vector* v,GIP_Vector* w,GIP_Vector* dxs,GIP_Vector* temp,double gamma,double rho)
{
    int threads=128;
    int blocks= (myTopo->getInnerSize() +threads -1)/threads;
    double ldt,dt;
    ldt=GIP_gpuGetCFL(u->getDevicePtr(),v->getDevicePtr(),w->getDevicePtr(),dxs->getDevicePtr(),temp->getHostPtr(),temp->getDevicePtr(),myTopo->getInnerSize(),gamma,rho,threads,blocks);
    
    MPI_Allreduce(&ldt,&dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
     return dt;

}

