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
        //GIP_Matrix(double* csrValA, int* csrColIndA, csrRowIndA,rows,cols,nnz) Link TF
        GIP_Matrix(char*&,GIP_Arch*);  // Read from a File in CSR
        virtual ~GIP_Matrix(); 
        virtual void postConstruct(char* &,GIP_Topo*,GIP_Arch*){};


        virtual void spmv(GIP_Vector*,GIP_Vector*);         // SpMV basic    y=Ax  
        virtual void spmv(GIP_Vector*,GIP_Vector*,double,double){};         // y=alpha*Ax+beta*y     
        virtual void spmv(GIP_Vector*,GIP_Matrix*){};         // SpMV that the result is stored in the csrVal of other matrix  y=Ax   
        virtual void spmv(GIP_Vector*,GIP_Matrix*,double,double){}; // SpMV that the result is stored in the csrVal of other matrix  y=alpha*Ax+beta*y         
        virtual void setUpDiagonal(GIP_Matrix*){};
        virtual void printRows(int); //Print some amount of Rows    
        virtual void operator=(GIP_Matrix){};


        virtual double* getCsrValADevice(){return NULL;};
        int getNumNnz(){return nnz;}; //temporal solo para pruebas
        int getNumRows(){return num_rows;}; //temporal solo para pruebas
        int getNumCols(){return num_cols;}; //temporal solo para pruebas
        double* getCsrValA(){return csrValA;}; //temporal solo para pruebas
        int* getCsrColIndA(){return csrColIndA;}; //temporal solo para pruebas
        int* getCsrRowIndA(){return csrRowIndA;}; //temporal solo para pruebas
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
