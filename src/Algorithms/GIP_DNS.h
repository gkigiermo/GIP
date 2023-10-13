#ifndef _gip_dns_
#define _gip_dns_

#include<stdio.h>
#include<iostream>
#include<string>
#include<mpi.h>
#include "GIP_LinearSolver.h"
#include "../BasicInterfaces/GIP_Topo.h"

using namespace std;

template<class Matrix,class Vector,class Node,class NonLinear> class GIP_DNS{
    public:
        GIP_DNS(){};
        GIP_DNS(char* name);
        ~GIP_DNS(){};

        void setUp();
        void integrate();

    private:

        string name;
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
    
        NonLinear noLin;

        // Parameters
        double Re;
        double gamma;
        double sigma;
        double rho;
        double dt;
        double time;
        int    maxIt;
        int    it;

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

};

template<class Matrix,class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::integrate()
{

   dt=noLin.getCFL(&topoCols,&u,&v,&w,&dxs,&temp,gamma,rho);
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

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::setUp() 
{
    uploadTopos();
    uploadConvDiff();
    uploadPoisson();
    uploadGradients();
    uploadMassCorrection();
    uploadCFL();

}

template<class Matrix,class Vector,class Node,class NonLinear> 
GIP_DNS<Matrix,Vector,Node,NonLinear>::GIP_DNS(char* _name):nodo(1,4,4),noLin(),topoCols(),topoFaces(),Ec(),D(),C(),L(),Lb(),
    Mx(),My(),Mz(),Gx(),Gy(),Gz(),SFx(),SFy(),SFz(),SO(),Mres(),dxs(),u(),v(),w(),u0(),v0(),w0(),ax(),ay(),
    az(),ax0(),ay0(),az0(),rhs(),p(),p0(),gpx(),gpy(),gpz(),cg(),temp(),temp2(){

        maxIt=10;
        rho=1.0;
        sigma=1.0/dt;
        dt=1e-3;
        Re=1.0e3;
        gamma=1.0/Re;
        time=0.0;
        name = _name;
    }

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::uploadConvDiff() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    string Dname = "InputFiles/Operators/CD/"+ name + "/" + name + "_" + to_string(nz) + "p/D"+name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    D.postConstruct(Dname,&topoCols,&nodo);
    C.postConstruct(Dname,&topoCols,&nodo);

    string Ecname = "InputFiles/Operators/CD/"+ name + "/" + name + "_" + to_string(nz) + "p/Ec" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Ec.postConstruct(Ecname,&topoFaces,&nodo);

}

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::uploadPoisson() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    string Lname = "InputFiles/Operators/L/"+ name + "/" + name + "_" + to_string(nz) + "p/L" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    L.postConstruct(Lname,&topoCols,&nodo);
    
    string Lbname = "InputFiles/Operators/L/"+ name + "/" + name + "_" + to_string(nz) + "p/Lb" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Lb.postConstruct(Lbname,&topoCols,&nodo);

    string Mxname = "InputFiles/Operators/M/"+ name + "/" + name + "_" + to_string(nz) + "p/Mx" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Mx.postConstruct(Mxname,&topoCols,&nodo);

    string Myname = "InputFiles/Operators/M/"+ name + "/" + name + "_" + to_string(nz) + "p/My" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    My.postConstruct(Myname,&topoCols,&nodo);

    string Mzname = "InputFiles/Operators/M/"+ name + "/" + name + "_" + to_string(nz) + "p/Mz" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Mz.postConstruct(Mzname,&topoCols,&nodo);

    GIP_Parameters param;
    param.setDouble("tol",1e-10);
    param.setBool("precond",true);
    param.setInt("maxIt",100);

    cg.setUp(&L,&param);

}

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::uploadGradients() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    string Gxname = "InputFiles/Operators/G/"+ name + "/" + name + "_" + to_string(nz) + "p/Gx" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Gx.postConstruct(Gxname,&topoCols,&nodo);

    string Gyname = "InputFiles/Operators/G/"+ name + "/" + name + "_" + to_string(nz) + "p/Gy" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Gy.postConstruct(Gyname,&topoCols,&nodo);

    string Gzname = "InputFiles/Operators/G/"+ name + "/" + name + "_" + to_string(nz) + "p/Gz" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Gz.postConstruct(Gzname,&topoCols,&nodo);

}

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::uploadMassCorrection() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    string SFxname = "InputFiles/Operators/S/"+ name + "/" + name + "_" + to_string(nz) + "p/SFx" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    SFx.postConstruct(SFxname,&topoFaces,&nodo);

    string SFyname = "InputFiles/Operators/S/"+ name + "/" + name + "_" + to_string(nz) + "p/SFy" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    SFy.postConstruct(SFyname,&topoFaces,&nodo);

    string SFzname = "InputFiles/Operators/S/"+ name + "/" + name + "_" + to_string(nz) + "p/SFz" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    SFz.postConstruct(SFzname,&topoFaces,&nodo);

    string SOname = "InputFiles/Operators/S/"+ name + "/" + name + "_" + to_string(nz) + "p/SO" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    SO.postConstruct(SOname,&topoFaces,&nodo);

}

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::uploadCFL()  
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    string Mresname = "InputFiles/Operators/M/"+ name + "/" + name + "_" + to_string(nz) + "p/Mres" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".csr"; 
    Mres.postConstruct(Mresname,&topoCols,&nodo);

}

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::uploadTopos() 
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    char Vecname[100];
    string Tname = "InputFiles/Topos/"+ name + "/" + name + "_" + to_string(nz) + "p/topoCols" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".topo"; 
    topoCols.postConstruct(Tname);

    string Fname = "InputFiles/Topos/"+ name + "/" + name + "_" + to_string(nz) + "p/topoFaces" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".topo"; 
    topoFaces.postConstruct(Fname);

    nodo.setUp();

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

    string VecUname = "InputFiles/Scalarfields/"+ name + "/" + name + "_" + to_string(nz) + "p/u" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".vec"; 
    u.postConstruct(VecUname,&topoCols,&nodo);
    u0.postConstruct(VecUname,&topoCols,&nodo);

    string VecVname = "InputFiles/Scalarfields/"+ name + "/" + name + "_" + to_string(nz) + "p/v" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".vec"; 
    v.postConstruct(VecVname,&topoCols,&nodo);
    v0.postConstruct(VecVname,&topoCols,&nodo);

    string VecWname = "InputFiles/Scalarfields/"+ name + "/" + name + "_" + to_string(nz) + "p/w" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".vec"; 
    w.postConstruct(VecWname,&topoCols,&nodo);
    w0.postConstruct(VecWname,&topoCols,&nodo);

    string VecDxname = "InputFiles/Scalarfields/"+ name + "/" + name + "_" + to_string(nz) + "p/dxs" + name + "_" + to_string(nz) + "p-" + to_string(rank) + ".vec"; 
    dxs.postConstruct(VecDxname,&topoCols,&nodo);

    ax=0;
    ay=0;
    az=0;
    ax0=0;
    ay0=0;
    az0=0;
    gpx=0;
    gpy=0;
    gpz=0;
    p=0;
    p0=0;
    mf=0; 
    rhs=0;
    temp=0;
    temp2=0;
   
    cout.setf(ios::scientific,ios::floatfield);

   }

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::PredictorVelocity()
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

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::PoissonEquation() 
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

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::VelocityCorrection() 
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

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::MassFluxes() 
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

template<class Matrix, class Vector,class Node,class NonLinear>
void GIP_DNS<Matrix,Vector,Node,NonLinear>::CFLCondition() 
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    Mres.spmv(&mf,&temp);
    double massres=temp.norm(&temp2,_INNER_);
  
    time+=dt; 
    if(rank==0) 
      cout<<"time "<<time<<" [dt::"<<dt<<"] massr "<<massres<<" ite "<<it+1<<endl;
    
    dt=noLin.getCFL(&topoCols,&u,&v,&w,&dxs,&temp,gamma,rho);

    sigma=1.0/dt;

    u.copyTo(&u0,_INNER_);
    v.copyTo(&v0,_INNER_);
    w.copyTo(&w0,_INNER_);

}


#endif
