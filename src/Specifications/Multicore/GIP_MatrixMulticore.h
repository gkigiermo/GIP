#ifndef _gip_matrixmulticore_
#define _gip_matrixmulticore_
#include <omp.h>
#include "../../AbstractInterfaces/GIP_Matrix.h"

class GIP_MatrixMulticore : public GIP_Matrix{
    public:
        GIP_MatrixMulticore(){};
        GIP_MatrixMulticore(char* name,GIP_Arch*);
        void postConstruct(string name,GIP_Topo*,GIP_Arch*);
        
        void spmv(GIP_Vector*,GIP_Vector*);
        void spmv(GIP_Vector*,GIP_Vector*,double,double);
        void spmv(GIP_Vector*,GIP_Matrix*);
        void spmv(GIP_Vector*,GIP_Matrix*,double,double);
        void setUpDiagonal(GIP_Matrix*);
        void printRows(int);
        ~GIP_MatrixMulticore();
        void operator=(GIP_Matrix&);
};

#endif
