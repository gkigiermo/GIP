#ifndef _gip_archcpu_
#define _gip_archcpu_

#define SetBitListOnMacro(A, B) { int pos, lpos; \
    pos = B / (sizeof(unsigned char) * 8); \
    lpos = B % (sizeof(unsigned char) * 8); \
    if( (A[pos] >> lpos) % 2 == 0 ) A[pos] += (1 << lpos); }

#define CheckBitListOnMacro(A, B, C) { int pos, lpos; \
    pos = B / (sizeof(unsigned char) * 8); \
    lpos = B % (sizeof(unsigned char) * 8); \
    C = (A[pos] >> lpos) % 2; }

#define SetBitListOffMacro(A, B) { int pos, lpos; \
    pos = B / (sizeof(unsigned char) * 8); \
    lpos = B % (sizeof(unsigned char) * 8); \
    if( (A[pos] >> lpos) % 2 == 1 ) A[pos] -= (1 << lpos); }

#include "../../AbstractInterfaces/GIP_Arch.h"
#include <mpi.h>
#include <omp.h>

class GIP_ArchCPU: public GIP_Arch{

    public:
        GIP_ArchCPU(){};
        GIP_ArchCPU(int, int , int);
        ~GIP_ArchCPU(){};

        void postConstruct(int,int,int);
        void setUp();

    protected:
        int devId;
        
};

#endif
