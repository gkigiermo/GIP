#ifndef _gip_vectorcuda_
#define _gip_vectorcuda_
#include "../../AbstractInterfaces/GIP_Vector.h"
#include "Headers/GIP_gpu.h"

class GIP_VectorCuda : public GIP_Vector{
    public:
        GIP_VectorCuda(){};
        GIP_VectorCuda(int,GIP_Arch*);
        GIP_VectorCuda(GIP_Topo*,GIP_Arch*);
        ~GIP_VectorCuda();

        void postConstruct(GIP_Topo*,GIP_Arch*);
        void postConstruct(char *,GIP_Topo*,GIP_Arch*);


        void axpy(GIP_Vector*, double,enum RUN_DOMAIN);
        void axpy(GIP_Vector*, double,double, enum RUN_DOMAIN);
        double dot(GIP_Vector*,GIP_Vector*, enum RUN_DOMAIN);
        double norm(GIP_Vector*, enum RUN_DOMAIN);
        void copyTo(GIP_Vector*, enum RUN_DOMAIN);
        void update();
        double* getDevicePtr();
        void FillRandom();
        void TransferToDevice(enum RUN_DOMAIN);
        void TransferToHost(enum RUN_DOMAIN);

        double* getHostPtr();
        void operator=(GIP_VectorCuda&);
        void operator=(int);
    protected:

        int threads;
        int blocks;
        double *dvec;
};
#endif
