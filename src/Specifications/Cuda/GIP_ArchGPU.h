#ifndef _gip_archgpu_
#define _gip_archgpu_
#include "../../AbstractInterfaces/GIP_Arch.h"
#include "Headers/GIP_gpu.h"
#include <mpi.h>

class GIP_ArchGPU: public GIP_Arch{

    public:
        GIP_ArchGPU(){};
        GIP_ArchGPU(int, int , int);
        ~GIP_ArchGPU(){};

        void postConstruct(int,int,int);
        void setUp();        

    protected:
        int devId;
};

#endif
