#include "GIP_VectorCuda.h"

GIP_VectorCuda::GIP_VectorCuda(int n,GIP_Arch* _node): GIP_Vector(n,_node) {

    GIP_cudaMallocDouble(dvec,n);
    GIP_cudaMemcpyDToGpu(vec,dvec,n);
    threads=128;
    blocks=(n+(threads-1))/threads;
}

GIP_VectorCuda::GIP_VectorCuda(GIP_Topo* topo,GIP_Arch* _node): GIP_Vector(topo,_node) {
    GIP_cudaMallocDouble(dvec,size);
    GIP_cudaMemcpyDToGpu(vec,dvec,size);
    threads=128;
    blocks=(size+(threads-1))/threads;
}
GIP_VectorCuda::~GIP_VectorCuda()
{
    if(dvec!=NULL)
        GIP_cudaFree(dvec);
}

void GIP_VectorCuda::postConstruct(GIP_Topo* topo,GIP_Arch* _node)
{
    size=topo->getAllSize();
    vec=new double[size];
    GIP_cudaMallocDouble(dvec,size);
    threads=128;
    blocks=(size+(threads-1))/threads;
    myTopo=topo;
    myArch=_node;
}

void GIP_VectorCuda::postConstruct(char * sfd_name,GIP_Topo* topo,GIP_Arch* _node)
{
    size=topo->getAllSize();
    vec=new double[size];
    GIP_cudaMallocDouble(dvec,size);
    threads=128;
    blocks=(size+(threads-1))/threads;
    myTopo=topo;
    myArch=_node;
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);
    FILE *fp;

    fp= fopen(sfd_name,"rb");

    int sizes[1];
/*
    for(int i=0;i<1;i++)
        fscanf(fp," %d",&sizes[0]);
    fscanf(fp," \n");

    for(int i=0;i<size;i++)
        fscanf(fp," %lf",&vec[i]);
 */
    fread(&sizes,sizeof(int),1,fp);
    fread(&vec[0],sizeof(double),sizes[0],fp);
	
    fclose(fp);
 
    GIP_cudaMemcpyDToGpu(vec,dvec,myTopo->getAllSize());
                 
}

void GIP_VectorCuda::TransferToDevice(enum RUN_DOMAIN rdom)
{
    switch(rdom)
    {
        case _ALL_:
                 GIP_cudaMemcpyDToGpu(vec,dvec,size);
                 break;
        case _OWNED_:
                 GIP_cudaMemcpyDToGpu(vec,dvec,myTopo->getOwnedSize());
                 break;
        case _INNER_:
                 GIP_cudaMemcpyDToGpu(vec,dvec,myTopo->getInnerSize());
                 break;
        case _HALO_:
                 GIP_cudaMemcpyDToGpu(vec+myTopo->getOwnedSize(),dvec+myTopo->getOwnedSize(),myTopo->getHaloSize());
                 break;
 
    };
}
void GIP_VectorCuda::TransferToHost(enum RUN_DOMAIN rdom)
{
     switch(rdom)
    {
        case _ALL_:
                  GIP_cudaMemcpyDToCpu(vec,dvec,size);
                     break;
       case _OWNED_:
                  GIP_cudaMemcpyDToCpu(vec,dvec,myTopo->getOwnedSize());
                     break;
       case _INNER_:
                  GIP_cudaMemcpyDToCpu(vec,dvec,myTopo->getInnerSize());
                     break;
       case _IFACE_:
                  GIP_cudaMemcpyDToCpu(vec+myTopo->getInterSize(),dvec+myTopo->getInterSize(),myTopo->getIfaceSize());
                     break;
    };
}
double* GIP_VectorCuda::getDevicePtr()
{
    return dvec;
}

void GIP_VectorCuda::axpy(GIP_Vector* x,double alpha,enum RUN_DOMAIN rdom)
{
    switch(rdom)
    {
        case _ALL_:
            GIP_cudaDaxpy(size,dvec,x->getDevicePtr(),alpha,blocks,threads);
            break;
        case _OWNED_:
            GIP_cudaDaxpy(myTopo->getOwnedSize(),dvec,x->getDevicePtr(),alpha,blocks,threads);
            break;
        case _INNER_:
            GIP_cudaDaxpy(myTopo->getInnerSize(),dvec,x->getDevicePtr(),alpha,blocks,threads);
            break; 
    };
}

void GIP_VectorCuda::axpy(GIP_Vector* x,double alpha,double beta,enum RUN_DOMAIN rdom)
{
    switch(rdom)
    {
        case _ALL_:
            GIP_cudaDaxpy(size,dvec,x->getDevicePtr(),alpha,beta,blocks,threads);
            break;
        case _OWNED_:
            GIP_cudaDaxpy(myTopo->getOwnedSize(),dvec,x->getDevicePtr(),alpha,beta,blocks,threads);
            break;
        case _INNER_:
            GIP_cudaDaxpy(myTopo->getInnerSize(),dvec,x->getDevicePtr(),alpha,beta,blocks,threads);
            break;
    };
}

double GIP_VectorCuda::dot(GIP_Vector* x, GIP_Vector* temp,enum RUN_DOMAIN rdom)
{
    int sz=MPI_Comm_size(MPI_COMM_WORLD,&sz);
    double tdot=0.0,dot=0.0;
    switch(rdom)
    {
        case _ALL_:
            tdot= GIP_cudaDdot(size,dvec,x->getDevicePtr(),temp->getHostPtr(),temp->getDevicePtr(),myTopo->getBlocksAll(),myTopo->getThreads());
            break;
        case _OWNED_:
            tdot= GIP_cudaDdot(myTopo->getOwnedSize(),dvec,x->getDevicePtr(),temp->getHostPtr(),temp->getDevicePtr(),myTopo->getBlocksOwned(),myTopo->getThreads());
            break;
        case _INNER_:
            tdot= GIP_cudaDdot(myTopo->getInnerSize(),dvec,x->getDevicePtr(),temp->getHostPtr(),temp->getDevicePtr(),myTopo->getBlocksInner(),myTopo->getThreads());
            break;
    };

    if(sz==1)
    {
        return tdot;
    }
    else
    {
        MPI_Allreduce(&tdot,&dot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        return dot;
    }
}

double GIP_VectorCuda::norm(GIP_Vector* temp,enum RUN_DOMAIN rdom)
{
    int sz=MPI_Comm_size(MPI_COMM_WORLD,&sz);
    double tnorm=0.0,norm=0.0;
    switch(rdom)
    {
        case _ALL_:
            tnorm= GIP_gpuNormi(dvec,temp->getHostPtr(),temp->getDevicePtr(),size,myTopo->getBlocksAll(),myTopo->getThreads());
            break;
        case _OWNED_:
            tnorm= GIP_gpuNormi(dvec,temp->getHostPtr(),temp->getDevicePtr(),myTopo->getOwnedSize(),myTopo->getBlocksOwned(),myTopo->getThreads());
 //           tdot= GIP_gpuNormi(myTopo->getOwnedSize(),dvec,x->getDevicePtr(),temp->getHostPtr(),temp->getDevicePtr(),myTopo->getBlocksOwned(),myTopo->getThreads());
            break;
        case _INNER_:
            tnorm= GIP_gpuNormi(dvec,temp->getHostPtr(),temp->getDevicePtr(),myTopo->getInnerSize(),myTopo->getBlocksInner(),myTopo->getThreads());
 //           tdot= GIP_gpuNormi(myTopo->getInnerSize(),dvec,x->getDevicePtr(),temp->getHostPtr(),temp->getDevicePtr(),myTopo->getBlocksInner(),myTopo->getThreads());
            break;
    };

    if(sz==1)
    {
        return tnorm;
    }
    else
    {
        MPI_Allreduce(&tnorm,&norm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        return norm;
    }
}
void GIP_VectorCuda::copyTo(GIP_Vector* x,enum RUN_DOMAIN rdom)
{
    switch(rdom)
    {
        case _ALL_:
             GIP_cudaDcopy(size,x->getDevicePtr(),dvec);
             break;
        case _OWNED_:
             GIP_cudaDcopy(myTopo->getOwnedSize(),x->getDevicePtr(),dvec);
             break;
        case _INNER_:
             GIP_cudaDcopy(myTopo->getInnerSize(),x->getDevicePtr(),dvec);
             break;
    };
}
void GIP_VectorCuda::operator=(GIP_VectorCuda& vector)
{
    size=vector.size;
    vec= new double[size];

    for(int i=0;i<size;i++)
        vec[i]=vector.vec[i];

    GIP_cudaMallocDouble(dvec,size);
    GIP_cudaMemcpyDToGpu(vec,dvec,size);
}

void GIP_VectorCuda::operator=(int value)
{
    for(int i=0;i<size;i++)
        vec[i]=value;
    GIP_cudaMemcpyDToGpu(vec,dvec,size);
}

void GIP_VectorCuda::update()
{
     double *send_buffer,*recv_buffer;
     send_buffer= new double[myTopo->getSendSize()];
     recv_buffer= new double[myTopo->getRecvSize()];
     int *send_indx=myTopo->getSendIndx();
     int *recv_indx=myTopo->getRecvIndx();

     GIP_cudaMemcpyDToCpu(vec+myTopo->getInterSize(),dvec+myTopo->getInterSize(),myTopo->getIfaceSize());
     
     for(int i=0;i<myTopo->getSendSize();i++)
     {
         send_buffer[i]=vec[send_indx[i]];
     }

     myTopo->sendrecv(send_buffer,recv_buffer);

     for(int i=0;i<myTopo->getRecvSize();i++)
     {
       vec[recv_indx[i]]=recv_buffer[i];
     }

     GIP_cudaMemcpyDToGpu(vec+myTopo->getOwnedSize(),dvec+myTopo->getOwnedSize(),myTopo->getHaloSize());
}
void GIP_VectorCuda::FillRandom()
{
    for(int i=0;i<size;i++)
       vec[i]=cos(i)+cos(10*i)+cos(100*i)+cos(1000*i);

    GIP_cudaMemcpyDToGpu(vec,dvec,size);
}

double* GIP_VectorCuda::getHostPtr()
{
//    GIP_cudaMemcpyDToCpu(vec,dvec,size);
    return vec;
}

