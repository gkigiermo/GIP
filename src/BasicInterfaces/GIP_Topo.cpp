#include "GIP_Topo.h"

void GIP_Topo::postConstruct(string filename)
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);

    FILE *fp;
    fp= fopen(filename.c_str(),"r");

    int sizes[9];
    for(int i=0;i<9;i++)
        fscanf(fp," %d",&sizes[i]);
    fscanf(fp," \n");

    send_size=sizes[0];    
    recv_size=sizes[1];
    cols_size=sizes[2];
    ownedsize=sizes[3];
    allsize=sizes[4];
    interiors=sizes[5];
    interfaces=sizes[6];
    boundaries=sizes[7];
    halos=sizes[8];
    innersize= interiors+interfaces;

    send_indx= new int[send_size];
    send_buffer= new double[send_size];
    recv_indx= new int[recv_size];
    recv_buffer= new double[recv_size];
    send_row= new int[cols_size];
    recv_row= new int[cols_size];

    for(int i=0;i<send_size;i++)
        fscanf(fp," %d",&send_indx[i]);
    fscanf(fp," \n");

    for(int i=0;i<cols_size;i++)
        fscanf(fp," %d",&send_row[i]);
    fscanf(fp," \n");

    for(int i=0;i<recv_size;i++)
        fscanf(fp," %d",&recv_indx[i]);
    fscanf(fp," \n");

    for(int i=0;i<cols_size;i++)
        fscanf(fp," %d",&recv_row[i]);
    fscanf(fp," \n");

    fclose(fp);

    threads=128;
    blocks_all=(allsize+(threads-1))/threads;
    blocks_owned=(ownedsize+(threads-1))/threads;
    blocks_inner=(innersize+(threads-1))/threads;
}

int GIP_Topo::getAllSize()
{
    return allsize;
}
int GIP_Topo::getOwnedSize()
{
    return ownedsize;
}
int GIP_Topo::getInnerSize()
{
    return innersize;

}
int GIP_Topo::getInterSize()
{
    return interiors;

}
int GIP_Topo::getIfaceSize()
{
    return interfaces;

}
int GIP_Topo::getHaloSize()
{
    return halos;

}
int GIP_Topo::getBoundarySize()
{
    return boundaries;

}

int GIP_Topo::getThreads()
{
    return threads;
}
int GIP_Topo::getBlocksAll()
{
    return blocks_all;
}
int GIP_Topo::getBlocksOwned()
{
    return blocks_owned;
}
int GIP_Topo::getBlocksInner()
{
    return blocks_inner;
}
void GIP_Topo::pack(double* & distvect)
{
    for(int i=0;i<send_size;i++)
    {
        send_buffer[i]=distvect[send_indx[i]];
    }
}
void GIP_Topo::unpack(double* &distvect)
{
    for(int i=0;i<recv_size;i++)
    {
        distvect[recv_indx[i]]=recv_buffer[i];
    }
}

void GIP_Topo::sendrecv()
{

    MPI_Request *rqs = new MPI_Request[cols_size-1]; // request send
    MPI_Request *rqr = new MPI_Request[cols_size-1]; // request received 

    int dim_send,dim_recv;
    for(int i=0;i<cols_size-1;i++)
    {
        dim_send=send_row[i+1]-send_row[i];
        dim_recv=recv_row[i+1]-recv_row[i];
        if(dim_send>0) MPI_Isend(&send_buffer[send_row[i]],dim_send,MPI_DOUBLE,i,0,MPI_COMM_WORLD, &rqs[i]);
        if(dim_recv>0) MPI_Irecv(&recv_buffer[recv_row[i]],dim_recv,MPI_DOUBLE,i,0,MPI_COMM_WORLD, &rqr[i]);
    }

 
    MPI_Status stat[2];
    for (int i=0;i<cols_size-1;i++) {
        if (send_row[i+1]-send_row[i] ) MPI_Wait(& rqs[i], &stat[0]);
        if (recv_row[i+1]-recv_row[i] ) MPI_Wait(& rqr[i], &stat[1]);
    }

    delete[] rqs;
    delete[] rqr;
}

void GIP_Topo::sendrecv(double *send_buffer_t,double* recv_buffer_t)
{

    MPI_Request *rqs = new MPI_Request[cols_size-1]; // request send
    MPI_Request *rqr = new MPI_Request[cols_size-1]; // request received 

    int dim_send,dim_recv;
    for(int i=0;i<cols_size-1;i++)
    {
        dim_send=send_row[i+1]-send_row[i];
        dim_recv=recv_row[i+1]-recv_row[i];
        if(dim_send>0) MPI_Isend(&send_buffer_t[send_row[i]],dim_send,MPI_DOUBLE,i,0,MPI_COMM_WORLD, &rqs[i]);
        if(dim_recv>0) MPI_Irecv(&recv_buffer_t[recv_row[i]],dim_recv,MPI_DOUBLE,i,0,MPI_COMM_WORLD, &rqr[i]);
    }

 
    MPI_Status stat[2];
    for (int i=0;i<cols_size-1;i++) {
        if (send_row[i+1]-send_row[i] ) MPI_Wait(& rqs[i], &stat[0]);
        if (recv_row[i+1]-recv_row[i] ) MPI_Wait(& rqr[i], &stat[1]);
    }
    delete[] rqs;
    delete[] rqr;
}


void GIP_Topo::update(double* distvect)
{
    pack(distvect);
    sendrecv();
    unpack(distvect);
}
GIP_Topo::~GIP_Topo()
{
/*    if(send_row!=NULL)
        delete[] send_row;
    if(recv_row!=NULL)
        delete[] recv_row;
    if(recv_buffer!=NULL)
        delete[] recv_buffer;
    if(send_buffer!=NULL)
        delete[] send_buffer;
    if(send_indx!=NULL)
        delete[] send_indx;
    if(recv_indx!=NULL)
        delete[] recv_indx;
*/
}

void GIP_Topo::printTopoInfo()
{
    cout<<" Topo Info  "<<endl;
    cout<<" inners     :"<<innersize<<endl;
    cout<<" interiors  :"<<interiors<<endl;
    cout<<" interfaces :"<<interfaces<<endl;
    cout<<" boundaries :"<<boundaries<<endl;
    cout<<" halos      :"<<halos<<endl;
    cout<<" Ownedsize  :"<<ownedsize<<endl;
    cout<<" Allsize    :"<<allsize<<endl;
}
