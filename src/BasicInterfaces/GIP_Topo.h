#ifndef _gip_topo_
#define _gip_topo_

#include<stdio.h>
#include<iostream>
#include<mpi.h>

/* Order of the elements within the topology    
      _____________________________
     |       |          | INTERIOR |
     |       |  INNER   |__________|
     | OWNED |          |INTERFACE |
     |       |__________|__________|
ALL->|       |BOUNDARIES|
     |_______|__________|
     | HALO  |   HALO   |
     |_______|__________|
*/
enum RUN_DOMAIN{
    _ALL_,   //TODOS
    _OWNED_, //OWNED 
    _HALO_,  //HALO
    _INNER_, // OWNED NOT BOUNDARY
    _INTER_, // INNER NOT IFACE
    _IFACE_, // INNER IFACE
    _BOUND_  // BOUNDARIES
};


using namespace std;

class GIP_Topo{

    protected:

        int allsize;
        int ownedsize;
        int innersize;

        int interiors;
        int interfaces;
        int boundaries;
        int halos;

        // Para comm;
        int send_size;
        int recv_size;
        int cols_size;

        int *send_indx;
        int *send_row;

        int *recv_indx;
        int *recv_row;


        ///////////////////////////////////
        //Estos variables deben ser redistribuidas en el vector
        double *send_buffer;
        double *recv_buffer;
        void pack(double *&);
        void unpack(double *&);
        void sendrecv(); // Prontamente tendra que ser eliminado para no basarse en ningun supuesto
        //////////////////////////////////
        

        int threads;
        int blocks_owned;
        int blocks_all;
        int blocks_inner;

    public:
        GIP_Topo(){};
        ~GIP_Topo();
        void postConstruct(char*);
        void update(double*); 
        int  getAllSize();
        int getOwnedSize();
        int getInnerSize();
        int getInterSize();
        int getIfaceSize();
        int getHaloSize();
        int getBoundarySize();

        int getSendSize(){return send_size;};
        int getRecvSize(){return recv_size;};
        int * getSendIndx(){return send_indx;};
        int * getSendRow(){return send_row;};
        int * getRecvIndx(){return recv_indx;};
        int * getRecvRow(){return recv_row;};

        void sendrecv(double*,double*); //Este es el procedimiento final



        // Esto deberia ir a una clase con info de la arquitectura
        int getThreads();
        int getBlocksAll();
        int getBlocksOwned();
        int getBlocksInner();
        ////////////////////////////////////////

        void printTopoInfo();
};

#endif
