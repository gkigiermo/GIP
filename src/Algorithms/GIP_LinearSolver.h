#ifndef _gip_linsolver_
#define _gip_linsolver_
#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<iostream>
#include "../BasicInterfaces/GIP_Topo.h"
#include "../BasicInterfaces/GIP_Parameters.h"
#include "../AbstractInterfaces/GIP_Arch.h"

template<class Matrix,class Vector> class GIP_LinearSolver{

    public:
        GIP_LinearSolver(){};
        ~GIP_LinearSolver(){};
        void setUp(Matrix*,GIP_Parameters*); 
        void solve(Vector*,Vector*,Vector*);
        double calcResid(Vector*,Vector*);
    protected: 

        Matrix* dA;
        Matrix diagonal;

        double tol;
        int maxIt;
        bool precond;
        GIP_Topo* myTopo;
        GIP_Arch* myNode;
};

template<class Matrix, class Vector>
void GIP_LinearSolver<Matrix,Vector>::setUp(Matrix* _A,GIP_Parameters* param)
{
    tol= param->getDouble("tol");
    maxIt=param->getInt("maxIt");
    precond=param->getBool("precond");

    dA=_A;
    myTopo= _A->getMyTopo();
    diagonal.setUpDiagonal(_A);
    myNode=_A->getMyArch(); 
}

template<class Matrix, class Vector>
void GIP_LinearSolver<Matrix,Vector>::solve(Vector* db, Vector* dx, Vector* dx0)
{
    double resIni;
    int sz,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&sz);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    double alpha, beta, r0, r1,dn;

    Vector dp(myTopo,myNode);
    dx0->update();
    dx0->copyTo(&dp,_ALL_);

    Vector dAx(myTopo,myNode);
    Vector dr(myTopo,myNode);
    Vector ds(myTopo,myNode);
    Vector dtemp(myTopo,myNode);

    dtemp=0;

    db->copyTo(&dr,_INNER_);
    dx0->copyTo(&dAx,_ALL_);
    dx0->copyTo(dx,_ALL_);

    dA->spmv(dx,&dAx);
    dr.axpy(&dAx,-1.0,_INNER_); 

    if(precond)
    {
        diagonal.spmv(&dr,&dp);
    }
    else
    {
        dr.copyTo(&dp,_INNER_);
    }
    r1=dr.dot(&dp,&dtemp,_INNER_);


    resIni=sqrt(r1);
    int k=1;
    while(k<=maxIt)
    {
        dp.update();

        dA->spmv(&dp,&dAx);

        dn=dp.dot(&dAx,&dtemp,_INNER_);

        alpha=r1 / dn;

        dx->axpy(&dp,alpha,_INNER_);

        dr.axpy(&dAx,-1.0*alpha,_INNER_); 

        if(precond){
            diagonal.spmv(&dr,&ds);
        }
        else
        {
            dr.copyTo(&ds,_INNER_); 
        }

        r0=r1;

        r1= dr.dot(&ds,&dtemp,_INNER_);  

        beta=r1 / r0;

        dp.axpy(&ds,1.0,beta,_INNER_); 

        k++;
    }
}

template<class Matrix,class Vector>
double GIP_LinearSolver<Matrix,Vector>::calcResid(Vector* dx,Vector* db)
{

    double res=0.0;
    Vector db2(myTopo,myNode);
    Vector dtemp(myTopo,myNode);
    dtemp=0.0;
    dx->update();
    dA->spmv(dx,&db2);
    db2.axpy(db,1.0,-1.0,_OWNED_);
    res=db2.dot(&db2,&dtemp,_OWNED_);

    return(sqrt(res));
}

#endif
