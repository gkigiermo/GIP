#ifndef _gip_vectormulticore_
#define _gip_vectormulticore_
#include "../../AbstractInterfaces/GIP_Vector.h"
#include <omp.h>
class GIP_VectorMulticore : public GIP_Vector{
    public:
        GIP_VectorMulticore(){};
        GIP_VectorMulticore(int,GIP_Arch*);
        GIP_VectorMulticore(GIP_Topo*,GIP_Arch*);
        ~GIP_VectorMulticore();

        void postConstruct(GIP_Topo*,GIP_Arch*);
        void postConstruct(char *n,GIP_Topo*,GIP_Arch*);

        void TransferToDevice(enum RUN_DOMAIN){};
        void TransferToHost(enum RUN_DOMAIN){};
        void axpy(GIP_Vector*, double,enum RUN_DOMAIN);
        void axpy(GIP_Vector*, double,double, enum RUN_DOMAIN);
        double dot(GIP_Vector*,GIP_Vector*, enum RUN_DOMAIN);
        double norm(GIP_Vector*,enum RUN_DOMAIN);
        void copyTo(GIP_Vector*, enum RUN_DOMAIN);
        void  update();
        double* getDevicePtr(){};
        void FillRandom();
        void operator=(GIP_VectorMulticore&);
        void operator=(int);
    protected:

};
#endif
