#include "../AbstractInterfaces/GIP_Vector.h"

GIP_Vector::GIP_Vector(int n,GIP_Arch* _node)
{
    vec=new double[n];
    size=n;
    myArch=_node;
}
GIP_Vector::GIP_Vector(GIP_Topo* topo,GIP_Arch* _node)
{
    size=topo->getAllSize();
    vec=new double[size];
    myTopo=topo; 
    myArch=_node;
}
GIP_Vector::~GIP_Vector()
{
/*    if(vec!=NULL)
        delete vec;
*/
}

void GIP_Vector::FillRandom()
{
    for(int i=0;i<size;i++)
       vec[i]=cos(1.0*i)+cos(10.0*i)+cos(100.0*i)+cos(1000.0*i);
}
void GIP_Vector::PrintRows(int n)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for(int i=1492;i<1492+n;i++)
        cout<<rank<<i<<": "<<vec[i]<<endl;
}

void GIP_Vector::operator=(double value)
{
    for(int i=0;i<size;i++)
       vec[i]=value;
}

double* GIP_Vector::getHostPtr()
{
    return vec;
}

void GIP_Vector::printPlainFile(char * opname)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    char filename[100];
    sprintf(filename,"%s-%d.vec",opname,rank);
    FILE *fp;
    double val;
    fp= fopen(filename,"w");

    for(int i=0;i<1;i++)
        fprintf(fp," %d ",size);
    fprintf(fp," \n");

    for(int i=0;i<size;i++)
    {
        val=vec[i];
        fprintf(fp,"%d %e \n",i,val);
    }
    fprintf(fp," \n");

    fclose(fp);



}
