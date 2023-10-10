#include "GIP_Matrix.h"

GIP_Matrix::GIP_Matrix(char* &matrix_name,GIP_Arch* _architecture)
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);
    char filename[100];
    sprintf(filename,"Operators/%s_%dp-%d.csr",matrix_name,nz,rank);
    FILE *fp;

    fp= fopen(filename,"r");

    int sizes[5];

    for(int i=0;i<5;i++)
        fscanf(fp," %d",&sizes[i]);
    fscanf(fp," \n");

    nnz=sizes[0];
    num_cols=sizes[1];
    num_rows=sizes[2];
        
    csrValA=new double[nnz];
    csrColIndA=new int[nnz];
    csrRowIndA= new int[num_rows+1];


    for(int i=0;i<nnz;i++)
        fscanf(fp," %lf",&csrValA[i]);
    fscanf(fp," \n");
    for(int i=0;i<nnz;i++)
        fscanf(fp," %d",&csrColIndA[i]);
    fscanf(fp," \n");
    for(int i=0;i<num_rows+1;i++)
        fscanf(fp," %d",&csrRowIndA[i]);
    fscanf(fp," \n");

    fclose(fp);

    myArch=_architecture;
}

GIP_Matrix::~GIP_Matrix()
{
/*
    cout<<"Destructor Matrix padre "<<endl;
    getchar();
    if(csrValA!=NULL)
        delete csrValA;
 
    cout<<"Libero csrValA "<<endl;
    getchar();
    if(csrColIndA!=NULL)
        delete[] csrColIndA;
    cout<<"Libero csrColIndA "<<endl;
    getchar();
    if(csrRowIndA!=NULL)
        delete[] csrRowIndA;
    cout<<"Libero csrRowIndA "<<endl;
    getchar();
  */
}

void GIP_Matrix::printRows(int n)
{
    cout<<"Matrix Info Clase Base:"<<endl;
    cout<<" ROWS : "<<num_rows<<" COLS : "<<num_cols<<" NNZ : "<<nnz<<endl;
    cout<<"Coefficients: "<<endl;
    for(int i=0;i<n;i++)
    {   cout<<" "<<i<<" :"; 
        int jstart=csrRowIndA[i];
        int jend=csrRowIndA[i+1];
        for(int j=jstart;j<jend;j++)
            cout<<" "<<csrValA[j];
        cout<<endl;
    }
}

void GIP_Matrix::spmv(GIP_Vector* x, GIP_Vector* b)
{
    double sum;
    const double *pval;
    int i,j,jb,je;

    pval = csrValA;
    double *vec_in=x->getHostPtr();
    double *vec_out=b->getHostPtr();

    for (i=0;i<num_rows;i++) {
        sum = 0;
        jb = csrRowIndA[i];
        je = csrRowIndA[i+1];
        for (j=jb;j<je;j++)
            sum +=  ((double)(vec_in[csrColIndA[j]])) * (*pval++);
        vec_out[i] = (double)sum;

    }
}

GIP_Topo* GIP_Matrix::getMyTopo()
{
    return myTopo;
}

GIP_Arch* GIP_Matrix::getMyArch()
{
    return myArch;
}
