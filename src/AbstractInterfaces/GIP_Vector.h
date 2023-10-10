#ifndef _gip_vector_
#define _gip_vector_

#include<mpi.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include "../BasicInterfaces/GIP_Topo.h"
#include "../AbstractInterfaces/GIP_Arch.h"


using namespace std;

class GIP_Vector{

    public:
        GIP_Vector(){};
        GIP_Vector(int,GIP_Arch*);
        GIP_Vector(GIP_Topo*,GIP_Arch*);
        virtual ~GIP_Vector();
 
        virtual void postConstruct(GIP_Topo*,GIP_Arch*){};
        virtual void postConstruct(char* &,GIP_Topo*,GIP_Arch*){};
 
        //operadores
        virtual void operator=(double);
        void operator=(GIP_Vector&);
  
      
        virtual double* getDevicePtr(){}; 
        double* getHostPtr();
        int getSize(){return size;};                
        virtual void axpy(GIP_Vector*, double){};
        virtual void axpy(GIP_Vector*, double,double){};
        virtual double dot(GIP_Vector*,GIP_Vector*){};
        virtual void copyTo(GIP_Vector*){};
        virtual void update(){};

        virtual void FillRandom(); 
        void PrintRows(int);         
	    void printPlainFile(char* ); 
        virtual void TransferToDevice(enum RUN_DOMAIN){}; 
        virtual void TransferToHost(enum RUN_DOMAIN){}; 

        GIP_Topo* getMyTopo(){return myTopo;};
        GIP_Arch* getMyArch(){return myArch;};


    protected:
        double *vec;
        int size;
        GIP_Topo* myTopo;
        GIP_Arch* myArch;
};

#endif
