#ifndef _gip_nlinearcuda_
#define _gip_nlinearcuda_

//#include"GIP_VectorCuda.h"
#include "../../AbstractInterfaces/GIP_NonLinear.h"
#include "Headers/GIP_gpu.h"


class GIP_NonLinearCuda : public GIP_NonLinear{

    public:
        GIP_NonLinearCuda(){};
        ~GIP_NonLinearCuda(){};
        double getCFL(GIP_Topo*,GIP_Vector*,GIP_Vector*,GIP_Vector*,GIP_Vector*,GIP_Vector*,double,double);
    protected:
};

#endif
