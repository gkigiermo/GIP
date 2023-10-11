#ifndef _gip_nlinear_
#define _gip_nlinear_

#include "GIP_Vector.h"

class GIP_NonLinear{

    public:
        GIP_NonLinear(){};
        ~GIP_NonLinear(){};
        virtual double getCFL(GIP_Topo*,GIP_Vector*,GIP_Vector*,GIP_Vector*,GIP_Vector*,GIP_Vector*,double,double){return 0;};
    protected:
    int l;
};

#endif
