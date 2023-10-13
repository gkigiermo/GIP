#ifndef _gip_matrix_
#define _gip_matrix_
#include<mpi.h>
#include<stdlib.h>
#include<iostream>
#include"GIP_Vector.h"
#include"GIP_Arch.h"

using namespace std;

class GIP_Matrix {

    public:
        GIP_Matrix(){};      
        GIP_Matrix(char*&,GIP_Arch*);  // Read from a File in CSR
        virtual ~GIP_Matrix(); 
        virtual void postConstruct(string name, GIP_Topo*, GIP_Arch*){};


        virtual void spmv(GIP_Vector*,GIP_Vector*);         // SpMV basic    y=Ax  
        virtual void spmv(GIP_Vector*,GIP_Vector*,double,double){};         // y=alpha*Ax+beta*y     
        virtual void spmv(GIP_Vector*,GIP_Matrix*){};         // SpMV that the result is stored in the csrVal of other matrix  y=Ax   
        virtual void spmv(GIP_Vector*,GIP_Matrix*,double,double){}; // SpMV that the result is stored in the csrVal of other matrix  y=alpha*Ax+beta*y         
        virtual void setUpDiagonal(GIP_Matrix*){};
        virtual void printRows(int); //Print some amount of Rows    
        virtual void operator=(GIP_Matrix){};


        virtual double* getCsrValADevice(){return NULL;};
        int getNumNnz(){return nnz;}; 
        int getNumRows(){return num_rows;};
        int getNumCols(){return num_cols;}; 
        double* getCsrValA(){return csrValA;};
        int* getCsrColIndA(){return csrColIndA;}; 
        int* getCsrRowIndA(){return csrRowIndA;}; 
        GIP_Topo* getMyTopo();
        GIP_Arch* getMyArch();
        
     protected:

        double *csrValA;
        int    *csrColIndA;
        int    *csrRowIndA;    

        int num_rows;
        int num_cols;
        int nnz;        
        
        GIP_Topo* myTopo;        
        GIP_Arch* myArch;
};
#endif
