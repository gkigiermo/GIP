#ifndef _gip_matrixcuda_
#define _gip_matrixcuda_

#include "../../AbstractInterfaces/GIP_Matrix.h"
#include "Headers/GIP_gpu.h"

class GIP_MatrixCuda : public GIP_Matrix{
    public:
        GIP_MatrixCuda(){};
        GIP_MatrixCuda(char* n,GIP_Arch*);
        void postConstruct(char *n,GIP_Topo*,GIP_Arch*);

        void spmv(GIP_Vector*,GIP_Vector*);
        void spmv(GIP_Vector*,GIP_Vector*,double,double);
        void spmv(GIP_Vector*,GIP_Matrix*);
        void spmv(GIP_Vector*,GIP_Matrix*,double,double);
        void setUpDiagonal(GIP_Matrix*);
        void printRows(int);
        ~GIP_MatrixCuda();
        void operator=(GIP_Matrix&);
        int threads;
        int blocks;

        protected:

        double *dcsrValA;
        int    *dcsrColIndA;
        int    *dcsrRowIndA;    
        double* getCsrValADevice(){return dcsrValA;};
};

#endif
