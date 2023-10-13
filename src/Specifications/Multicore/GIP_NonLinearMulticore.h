#ifndef _gip_nlinearmulticore_
#define _gip_nlinearmulticore_

#include "../../AbstractInterfaces/GIP_NonLinear.h"


class GIP_NonLinearMulticore : public GIP_NonLinear{

    public:
        GIP_NonLinearMulticore(){};
        ~GIP_NonLinearMulticore(){};
        double getCFL(GIP_Topo*,GIP_Vector*,GIP_Vector*,GIP_Vector*,GIP_Vector*,GIP_Vector*,double,double);
    protected:
};

#endif
