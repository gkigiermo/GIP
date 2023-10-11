#include "GIP_VectorMulticore.h"

GIP_VectorMulticore::GIP_VectorMulticore(int n,GIP_Arch* _architecture): GIP_Vector(n,_architecture) {}

GIP_VectorMulticore::GIP_VectorMulticore(GIP_Topo* topo,GIP_Arch* _architecture): GIP_Vector(topo,_architecture) {}
GIP_VectorMulticore::~GIP_VectorMulticore(){}

void GIP_VectorMulticore::postConstruct(GIP_Topo* topo,GIP_Arch* _architecture)
{
    size=topo->getAllSize();
    vec=new double[size];
    myTopo=topo;
    myArch=_architecture;
}

void GIP_VectorMulticore::postConstruct(string sfd_name,GIP_Topo* topo,GIP_Arch* _architecture)
{

    size=topo->getAllSize();
    vec=new double[size];
    myTopo=topo;
    myArch=_architecture;

    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);
    FILE *fp;
    fp= fopen(sfd_name.c_str(),"rb");
    if (fp == NULL) {
    
        cout<<" Error opening scalarfield file "<<sfd_name<<endl;
        exit(0);
    }
    int tam=0;

    fread(&tam,sizeof(int),1,fp);
  
    fread(vec,sizeof(double),tam,fp);
    
    fclose(fp);
 
}

void GIP_VectorMulticore::axpy(GIP_Vector* x,double alpha,enum RUN_DOMAIN rdom)
{

    int axpy_size=0;
    double* xptr=x->getHostPtr();
    switch(rdom)
    {
        case _ALL_:
            axpy_size=size;
            break;
        case _OWNED_:
            axpy_size=myTopo->getOwnedSize();
            break;
        case _INNER_:
            axpy_size=myTopo->getInnerSize();
            break; 
    };

//    #pragma omp parallel for schedule(static)  
    for(int i=0;i<axpy_size;i++)
        vec[i]+=alpha*xptr[i];


}

void GIP_VectorMulticore::axpy(GIP_Vector* x,double alpha,double beta,enum RUN_DOMAIN rdom)
{
    int axpy_size=0;
    double* xptr=x->getHostPtr();
    switch(rdom)
    {
        case _ALL_:
            axpy_size=size;
            break;
        case _OWNED_:
            axpy_size=myTopo->getOwnedSize();
            break;
        case _INNER_:
            axpy_size=myTopo->getInnerSize();
            break; 
    };

//    #pragma omp parallel for schedule(static)  
    for(int i=0;i<axpy_size;i++)
        vec[i]=alpha*xptr[i]+beta*vec[i];



}
double GIP_VectorMulticore::dot(GIP_Vector* x, GIP_Vector* temp,enum RUN_DOMAIN rdom)
{
    temp = NULL; //UNUSED

    int sz=MPI_Comm_size(MPI_COMM_WORLD,&sz);
    double tdot=0.0,dot=0.0;

    int dot_size=0;
    double* xptr=x->getHostPtr();
    switch(rdom)
    {
        case _ALL_:
            dot_size=size;
            break;
        case _OWNED_:
            dot_size=myTopo->getOwnedSize();
            break;
        case _INNER_:
            dot_size=myTopo->getInnerSize();
            break; 
    };

//    #pragma omp parallel for reduction(+:tdot)
    for(int i=0;i<dot_size;i++)
            tdot+=vec[i]*xptr[i];

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

double GIP_VectorMulticore::norm(GIP_Vector* temp,enum RUN_DOMAIN rdom)
{
    temp = NULL; //UNUSED

    int sz=MPI_Comm_size(MPI_COMM_WORLD,&sz);
    double norm=0.0,tnorm=0.0;

    int norm_size=0;
    switch(rdom)
    {
        case _ALL_:
            norm_size=size;
            break;
        case _OWNED_:
            norm_size=myTopo->getOwnedSize();
            break;
        case _INNER_:
            norm_size=myTopo->getInnerSize();
            break; 
    };

    for(int i=0;i<norm_size;i++)
        if(fabs(vec[i])>tnorm)
            tnorm=fabs(vec[i]);

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


void GIP_VectorMulticore::copyTo(GIP_Vector* x,enum RUN_DOMAIN rdom)
{
    switch(rdom)
    {
        case _ALL_:
             memcpy(x->getHostPtr(),vec,sizeof(double)*size);
             break;
        case _OWNED_:
             memcpy(x->getHostPtr(),vec,sizeof(double)*myTopo->getOwnedSize());
             break;
        case _INNER_:
             memcpy(x->getHostPtr(),vec,sizeof(double)*myTopo->getInnerSize());
             break;
    };
}
void GIP_VectorMulticore::operator=(GIP_VectorMulticore& vector)
{
    size=vector.size;
    vec= new double[size];
 
    for(int i=0;i<size;i++)
        vec[i]=vector.vec[i];
}

void GIP_VectorMulticore::operator=(int value)
{
    for(int i=0;i<size;i++)
        vec[i]=value;
}

void GIP_VectorMulticore::update()
{
     double *send_buffer,*recv_buffer;
     send_buffer= new double[myTopo->getSendSize()];
     recv_buffer= new double[myTopo->getRecvSize()];
     int *send_indx=myTopo->getSendIndx();
     int *recv_indx=myTopo->getRecvIndx();
     
     for(int i=0;i<myTopo->getSendSize();i++)
     {
         send_buffer[i]=vec[send_indx[i]];
     }

     myTopo->sendrecv(send_buffer,recv_buffer);

     for(int i=0;i<myTopo->getRecvSize();i++)
     {
       vec[recv_indx[i]]=recv_buffer[i];
     }
}
void GIP_VectorMulticore::FillRandom()
{
    for(int i=0;i<size;i++)
       vec[i]=cos(1.0*i)+cos(10.0*i)+cos(100.0*i)+cos(1000.0*i);
}
