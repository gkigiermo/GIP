#ifndef _gip_dns_
#define _gip_dns_

#include<stdio.h>
#include<iostream>
#include<mpi.h>
#include "GIP_LinearSolver.h"
#include "../BasicInterfaces/GIP_Topo.h"
#include "../Specifications/Multicore/GIP_ArchCPU.h"


using namespace std;

template<class Matrix,class Vector,class Node> class GIP_DNS{


    public:

        GIP_DNS(){};
        GIP_DNS(char* name);
        ~GIP_DNS(){};

        void setUp();
        void integrate();

    private:

        char name[100]; 

        Node nodo;

        //Convection - Diffusion 
        Matrix C,Ec;
        Matrix D;

        //Poisson Equation
        Matrix L;
        Matrix Lb;  // To extend the boundary conditions 

        //Divergence and Gradient
        Matrix Mx,My,Mz;       
        Matrix Gx,Gy,Gz;

        //Kinematic energy correction 
        Matrix SFx,SFy,SFz,SO;

        //Mass flow operator
        Matrix Mres;

        //Scalarfields
        
        Vector u,v,w;      // Velocities step n
        Vector u0,v0,w0;   // Velocities step n-1
        
        Vector ax,ay,az;    // Adams Bashford n
        Vector ax0,ay0,az0; // Adams Bashford n-1

        Vector rhs;   // right hand size Poisson
        Vector p,p0;  // pressures n, n-1

        Vector gpx,gpy,gpz; //gradients of p in x,y,z

        Vector mf; // mass fluxes in the faces vector

        Vector dxs; // cfl deltaxs
        Vector temp;
        Vector temp2;


        // Linear Solver
        GIP_LinearSolver<Matrix,Vector> cg;

        // Topologies
        GIP_Topo topoCols;
        GIP_Topo topoFaces;



        // Parameters
        double Re;
        double gamma;
        double sigma;
        double rho;
        double dt;
        double time;
	int    maxIt;
        int    it;
	
        // Private methods
        void uploadTopos();  
        void uploadConvDiff(); 
        void uploadPoisson(); 
        void uploadGradients(); 
        void uploadMassCorrection(); 
        void uploadCFL(); 

        void PredictorVelocity(); 
        void PoissonEquation(); 
        void VelocityCorrection(); 
        void MassFluxes(); 
        void CFLCondition(); 

        void calculateCFL(); 
};

template<class Matrix,class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::integrate()
{

    calculateCFL();
    sigma=1.0/dt;
    for(it=0;it<maxIt;it++)
    {
        PredictorVelocity();
        PoissonEquation();
        VelocityCorrection();
        MassFluxes();
        CFLCondition();
    }

}


template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::setUp() 
{
    uploadTopos();
    uploadConvDiff();
    uploadPoisson();
    uploadGradients();
    uploadMassCorrection();
    uploadCFL();
    calculateCFL();
}

template<class Matrix,class Vector,class Node> 
GIP_DNS<Matrix,Vector,Node>::GIP_DNS(char* _name):nodo(1,4,4),topoCols(),topoFaces(),Ec(),D(),C(),L(),Lb(),
    Mx(),My(),Mz(),Gx(),Gy(),Gz(),SFx(),SFy(),SFz(),SO(),Mres(),dxs(),u(),v(),w(),u0(),v0(),w0(),ax(),ay(),
    az(),ax0(),ay0(),az0(),rhs(),p(),p0(),gpx(),gpy(),gpz(),cg(),temp(),temp2(){

        maxIt=100;
        rho=1.0;
        sigma=1.0/dt;
        dt=1e-3;
        Re=1.0e3;
        gamma=1.0/Re;
        time=0.0;
	sprintf(name,"%s",_name);
    }




template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::uploadConvDiff() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    char Dname[100];
    sprintf(Dname,"InputFiles/Operators/CD/%s/%s_%dp/D%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    D.postConstruct(&Dname[0],&topoCols,&nodo);
    C.postConstruct(&Dname[0],&topoCols,&nodo);

    char Ecname[100];
    sprintf(Ecname,"InputFiles/Operators/CD/%s/%s_%dp/Ec%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Ec.postConstruct(&Ecname[0],&topoFaces,&nodo);

}


template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::uploadPoisson() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    char Lname[100];
    sprintf(Lname,"InputFiles/Operators/L/%s/%s_%dp/L%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    L.postConstruct(&Lname[0],&topoCols,&nodo);

    char Lbname[100];
    sprintf(Lbname,"InputFiles/Operators/L/%s/%s_%dp/Lb%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Lb.postConstruct(&Lbname[0],&topoCols,&nodo);

    char Mxname[100];
    sprintf(Mxname,"InputFiles/Operators/M/%s/%s_%dp/Mx%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Mx.postConstruct(&Mxname[0],&topoCols,&nodo);

    char Myname[100];
    sprintf(Myname,"InputFiles/Operators/M/%s/%s_%dp/My%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    My.postConstruct(&Myname[0],&topoCols,&nodo);

    char Mzname[100];
    sprintf(Mzname,"InputFiles/Operators/M/%s/%s_%dp/Mz%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Mz.postConstruct(&Mzname[0],&topoCols,&nodo);


    GIP_Parameters param;
    param.setDouble("tol",1e-10);
    param.setBool("precond",true);
    param.setInt("maxIt",100);

    cg.setUp(&L,&param);

}


template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::uploadGradients() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);


    char Gxname[100];
    sprintf(Gxname,"InputFiles/Operators/G/%s/%s_%dp/Gx%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Gx.postConstruct(&Gxname[0],&topoCols,&nodo);

    char Gyname[100];
    sprintf(Gyname,"InputFiles/Operators/G/%s/%s_%dp/Gy%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Gy.postConstruct(&Gyname[0],&topoCols,&nodo);

    char Gzname[100];
    sprintf(Gzname,"InputFiles/Operators/G/%s/%s_%dp/Gz%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Gz.postConstruct(&Gzname[0],&topoCols,&nodo);

}


template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::uploadMassCorrection() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    char SFxname[100];
    sprintf(SFxname,"InputFiles/Operators/S/%s/%s_%dp/SFx%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    SFx.postConstruct(&SFxname[0],&topoFaces,&nodo);

    char SFyname[100];
    sprintf(SFyname,"InputFiles/Operators/S/%s/%s_%dp/SFy%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    SFy.postConstruct(&SFyname[0],&topoFaces,&nodo);

    char SFzname[100];
    sprintf(SFzname,"InputFiles/Operators/S/%s/%s_%dp/SFz%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    SFz.postConstruct(&SFzname[0],&topoFaces,&nodo);

    char SOname[100];
    sprintf(SOname,"InputFiles/Operators/S/%s/%s_%dp/SO%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    SO.postConstruct(&SOname[0],&topoFaces,&nodo);

}

template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::uploadCFL()  
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    char Mresname[100];
    sprintf(Mresname,"InputFiles/Operators/M/%s/%s_%dp/Mres%s_%dp-%d.csr",name,name,nz,name,nz,rank);
    Mres.postConstruct(&Mresname[0],&topoCols,&nodo);

}
template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::uploadTopos() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    char Tname[100];
    char Vecname[100];
    sprintf(Tname,"InputFiles/Topos/%s/%s_%dp/topoCols%s_%dp-%d.topo",name,name,nz,name,nz,rank);
    topoCols.postConstruct(&Tname[0]);
    sprintf(Tname,"InputFiles/Topos/%s/%s_%dp/topoFaces%s_%dp-%d.topo",name,name,nz,name,nz,rank);
    topoFaces.postConstruct(&Tname[0]);

    //nodo.setUp();

    p.postConstruct(&topoCols,&nodo);
    p0.postConstruct(&topoCols,&nodo);
    ax.postConstruct(&topoCols,&nodo);
    ax0.postConstruct(&topoCols,&nodo);
    ay.postConstruct(&topoCols,&nodo);
    ay0.postConstruct(&topoCols,&nodo);
    az.postConstruct(&topoCols,&nodo);
    az0.postConstruct(&topoCols,&nodo);
    gpx.postConstruct(&topoCols,&nodo);
    gpy.postConstruct(&topoCols,&nodo);
    gpz.postConstruct(&topoCols,&nodo);
    rhs.postConstruct(&topoCols,&nodo);
    dxs.postConstruct(&topoCols,&nodo);
    temp.postConstruct(&topoCols,&nodo);
    temp2.postConstruct(&topoCols,&nodo);
    mf.postConstruct(&topoFaces,&nodo);

    sprintf(Vecname,"InputFiles/Scalarfields/%s/%s_%dp/u%s_%dp-%d.vec",name,name,nz,name,nz,rank);
    u.postConstruct(Vecname,&topoCols,&nodo);
    u0.postConstruct(Vecname,&topoCols,&nodo);

    sprintf(Vecname,"InputFiles/Scalarfields/%s/%s_%dp/v%s_%dp-%d.vec",name,name,nz,name,nz,rank);
    v.postConstruct(Vecname,&topoCols,&nodo);
    v0.postConstruct(Vecname,&topoCols,&nodo);

    sprintf(Vecname,"InputFiles/Scalarfields/%s/%s_%dp/w%s_%dp-%d.vec",name,name,nz,name,nz,rank);
    w.postConstruct(Vecname,&topoCols,&nodo);
    w0.postConstruct(Vecname,&topoCols,&nodo);

    sprintf(Vecname,"InputFiles/Scalarfields/%s/%s_%dp/dxs%s_%dp-%d.vec",name,name,nz,name,nz,rank);
    dxs.postConstruct(Vecname,&topoCols,&nodo);

    ax=0.0;
    ay=0.0;
    az=0.0;
    ax0=0.0;
    ay0=0.0;
    az0=0.0;
    gpx=0.0;
    gpy=0.0;
    gpz=0.0;
    p=0.0;
    p0=0.0;
    mf=0.0; 
    rhs=0.0;
    temp=0.0;
    temp2=0.0;
   
    cout.setf(ios::scientific,ios::floatfield);


}


template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::PredictorVelocity()
{

    Ec.spmv(&mf,&C);  
    D.spmv(&u,&ax);
    D.spmv(&v,&ay);
    D.spmv(&w,&az);
    C.spmv(&u,&ax,-1.0,1.0);
    C.spmv(&v,&ay,-1.0,1.0);
    C.spmv(&w,&az,-1.0,1.0);

    double alpha1 = 1.5*dt;    
    double alpha2 = -0.5*dt;
    
    u.axpy(&ax,alpha1,_INNER_);
    u.axpy(&ax0,alpha2,_INNER_);
    v.axpy(&ay,alpha1,_INNER_);
    v.axpy(&ay0,alpha2,_INNER_);
    w.axpy(&az,alpha1,_INNER_);
    w.axpy(&az0,alpha2,_INNER_); 

    u.update();
    v.update();
    w.update();
    
    ax.copyTo(&ax0,_INNER_);
    ay.copyTo(&ay0,_INNER_);
    az.copyTo(&az0,_INNER_);

}
template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::PoissonEquation() 
{
    Mx.spmv(&u,&rhs);
    My.spmv(&v,&rhs,1.0,1.0);
    Mz.spmv(&w,&rhs,sigma,sigma);

    cg.solve(&rhs,&p,&p0);
    p.update();

    p.copyTo(&p0,_INNER_);
    p0.copyTo(&p,_OWNED_);
    Lb.spmv(&p0,&p,1.0,1.0);

}

template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::VelocityCorrection() 
{
    Gx.spmv(&p,&gpx);
    Gy.spmv(&p,&gpy);
    Gz.spmv(&p,&gpz);
    gpx.update();
    gpy.update();
    gpz.update();

    double fact= -1.0/(rho*sigma);
    
    u.axpy(&gpx,fact,_INNER_);
    v.axpy(&gpy,fact,_INNER_);
    w.axpy(&gpz,fact,_INNER_);

    u.update();
    v.update();
    w.update();
    
   
}

template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::MassFluxes() 
{
	double fact1=-1.0/sigma;
	double fact2=1.0/sigma;

	   SO.spmv(&p,&mf,fact1,0.0);
	   SFx.spmv(&gpx,&mf,fact2,1.0);
	   SFy.spmv(&gpy,&mf,fact2,1.0);
	   SFz.spmv(&gpz,&mf,fact2,1.0);
	    SFx.spmv(&u,&mf,rho,1.0);
	    SFy.spmv(&v,&mf,rho,1.0);
	    SFz.spmv(&w,&mf,rho,1.0);
}

template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::CFLCondition() 
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    Mres.spmv(&mf,&temp);
    double massdot=temp.dot(&temp,&temp,_OWNED_);
    double massres=temp.norm(&temp2,_OWNED_);
  
    time+=dt; 
    if(rank==0) 
      cout<<"time "<<time<<" [dt::"<<dt<<"] massr "<<massres<<" ite "<<it+1<<endl;
    
    calculateCFL();

    sigma=1.0/dt;

    u.copyTo(&u0,_INNER_);
    v.copyTo(&v0,_INNER_);
    w.copyTo(&w0,_INNER_);

}

template<class Matrix, class Vector,class Node>
void GIP_DNS<Matrix,Vector,Node>::calculateCFL() 
{

  //  dt=1e-3;
    dt=1e20;
    double convCoef=0.25;
    double difCoef=0.2;
    double ut,vt,wt,dx;
    double normit=0.0;
    double* uptr,*vptr,*wptr,*dxsptr;
    double aux;
    double limite=1e-17;
    uptr=u.getHostPtr();
    vptr=v.getHostPtr();
    wptr=w.getHostPtr();
    dxsptr=dxs.getHostPtr();

    for(int i=0;i<topoCols.getInnerSize();i++)
    {
        ut=uptr[i];
        vt=vptr[i];
        wt=wptr[i];
        dx=dxsptr[i];
        if(fabs(ut)>fabs(vt))
            normit=fabs(ut);
        else
            normit=fabs(vt);

        if(normit<fabs(wt))
            normit=fabs(wt);
        if(!(dx/normit < 1e-17))
        {
            if(dt>convCoef*dx/(normit))
                dt= convCoef*dx/(normit);
        }
        if(!(dx*dx*rho/gamma  < 1e-17))
        {
            if(dt>0.2*dx*dx*rho/gamma)
                dt=0.2*dx*dx*rho/gamma;
        }

    }
    double final = 0.0;
    MPI_Allreduce(&dt,&final,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);


    if(0.8*final>1e-1)
        dt=1e-1;
    else
        dt=0.8*final;

}
#endif
