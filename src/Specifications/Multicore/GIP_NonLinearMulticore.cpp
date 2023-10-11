#include "GIP_NonLinearMulticore.h"

double GIP_NonLinearMulticore::getCFL(GIP_Topo* myTopo,GIP_Vector* u,GIP_Vector* v,GIP_Vector* w,GIP_Vector* dxs,GIP_Vector* temp,double gamma,double rho)
{
    double dt;
   
    dt=1e20;
    double convCoef=0.25;
    double ut,vt,wt,dx;
    double normit=0.0;
    double* uptr,*vptr,*wptr,*tptr,*dxsptr;
    
    uptr=u->getHostPtr();
    vptr=v->getHostPtr();
    wptr=w->getHostPtr();
    tptr=temp->getHostPtr();
    dxsptr=dxs->getHostPtr();

    for(int i=0;i<myTopo->getInnerSize();i++)
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
   
    return dt;

}

