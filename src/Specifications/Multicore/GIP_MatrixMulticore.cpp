#include "GIP_MatrixMulticore.h"

GIP_MatrixMulticore::GIP_MatrixMulticore(char* n,GIP_Arch* _architecture): GIP_Matrix(n,_architecture) {}

void GIP_MatrixMulticore::postConstruct(char *matrix_name,GIP_Topo* _topo,GIP_Arch* _architecture)
{
    int rank,nz;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nz);
    FILE *fp;

    fp= fopen(matrix_name,"rb");
    int sizes[5];
    fread(&sizes,sizeof(int),5,fp);
    nnz=sizes[0];
    num_cols=sizes[1];
    num_rows=sizes[2];
    csrValA=new double[nnz];
    csrColIndA=new int[nnz];
    csrRowIndA= new int[num_rows+1];


    fread(csrValA,sizeof(double),nnz,fp);
    fread(csrColIndA,sizeof(int),nnz,fp);
    fread(csrRowIndA,sizeof(int),num_rows+1,fp);

    fclose(fp);
    /*
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
     */ 
    myTopo=_topo;
    myArch=_architecture;
}


GIP_MatrixMulticore::~GIP_MatrixMulticore()
{
    /*
       if(csrValA!=NULL)
       delete csrValA;
       if(csrColIndA!=NULL)
       delete[] csrColIndA;
       if(csrRowIndA!=NULL)
       delete[] csrRowIndA;
     */
}

void GIP_MatrixMulticore::printRows(int n)
{
    int nrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&nrank);
    cout<<"Matrix Info Clase CUDA :"<<endl;
    cout<<" ROWS : "<<num_rows<<" COLS : "<<num_cols<<" NNZ : "<<nnz<<endl;
    cout<<"Coefficients: "<<endl;
    for(int i=1404;i<1404+n;i++)
    {  int jstart=csrRowIndA[i];
        int jend=csrRowIndA[i+1];
         cout<<nrank<<": "<<i<<" :"<<jstart<<" "<<jend<<endl; 
        for(int j=jstart;j<jend;j++)
            cout<<" "<<csrValA[j]<<" "<<csrColIndA[j];
        cout<<endl;
    }
}


void GIP_MatrixMulticore::spmv(GIP_Vector* x, GIP_Vector* b)
{
    double sum;
    int i,j,jb,je;

    double *vec_in=x->getHostPtr();
    double *vec_out=b->getHostPtr();


//#pragma omp parallel for schedule(static) private(j,jb,je,sum)
    for (i=0;i<num_rows;i++) {
        sum = 0;
        jb = csrRowIndA[i];
        je = csrRowIndA[i+1];
        for (j=jb;j<je;j++)
            sum+=vec_in[csrColIndA[j]]*csrValA[j];
        if(jb!=je)
                 vec_out[i] = sum;
    }
}

void GIP_MatrixMulticore::spmv(GIP_Vector* x, GIP_Vector* b,double alpha,double beta)
{
    double sum;
    int i,j,jb,je;

    double *vec_in=x->getHostPtr();
    double *vec_out=b->getHostPtr();

//#pragma omp parallel for schedule(static) private(j,jb,je,sum)
    for (i=0;i<num_rows;i++) {
        sum = 0;
        jb = csrRowIndA[i];
        je = csrRowIndA[i+1];
        for (j=jb;j<je;j++)
            sum+=vec_in[csrColIndA[j]]*csrValA[j];
        if(jb!=je)    
          vec_out[i] = alpha*sum+beta*vec_out[i];
    }
}

void GIP_MatrixMulticore::spmv(GIP_Vector* x, GIP_Matrix* Ab)
{
    double sum;
    int i,j,jb,je;

    double *vec_in=x->getHostPtr();
    double *vec_out=Ab->getCsrValA();

//#pragma omp parallel for schedule(static) private(j,jb,je,sum)
    for (i=0;i<num_rows;i++) {
        sum = 0;
        jb = csrRowIndA[i];
        je = csrRowIndA[i+1];
        for (j=jb;j<je;j++)
            sum+=vec_in[csrColIndA[j]]*csrValA[j];
        if(jb!=je)
            vec_out[i] = sum;
    }
}

void GIP_MatrixMulticore::spmv(GIP_Vector* x, GIP_Matrix* Ab,double alpha,double beta)
{
    double sum;
    int i,j,jb,je;

    double *vec_in=x->getHostPtr();
    double *vec_out=Ab->getCsrValA();

//#pragma omp parallel for schedule(static) private(j,jb,je,sum)
    for (i=0;i<num_rows;i++) {
        sum = 0;
        jb = csrRowIndA[i];
        je = csrRowIndA[i+1];
        for (j=jb;j<je;j++)
            sum+=vec_in[csrColIndA[j]]*csrValA[j];
        if(jb!=je)
          vec_out[i] = alpha*sum+beta*vec_out[i];
    }
}

void GIP_MatrixMulticore::operator=(GIP_Matrix& matrix)
{
    num_rows=matrix.getNumRows();
    num_cols=matrix.getNumCols();
    nnz=matrix.getNumNnz();


    csrValA=new double[nnz];
    csrColIndA=new int[nnz];
    csrRowIndA= new int[num_rows+1];

    double* valA=matrix.getCsrValA();
    int*    colA=matrix.getCsrColIndA();
    int*    rowA=matrix.getCsrRowIndA();


    for(int i=0;i<nnz;i++)
    {
        csrValA[i]=valA[i];
        csrColIndA[i]=colA[i];
    }

    for(int i=0;i<num_rows+1;i++)
    {
        csrRowIndA[i]=rowA[i];
    }

    myTopo=matrix.getMyTopo();
    myArch=matrix.getMyArch();

}

void GIP_MatrixMulticore::setUpDiagonal(GIP_Matrix* _A)
{
    num_rows=_A->getNumRows();
    num_cols=_A->getNumRows();
    nnz=_A->getNumRows();
    double* valA=_A->getCsrValA();
    int*    colA=_A->getCsrColIndA();
    int*    rowA=_A->getCsrRowIndA();

    csrValA=new double[nnz];
    csrColIndA=new int[nnz];
    csrRowIndA= new int[num_rows+1];


    for(int i=0;i<num_rows;i++)
    {
        int js=rowA[i];
        int je=rowA[i+1];
        for(int j=js;j<je;j++)
        {
            if(i==colA[j])
            {
                csrValA[i]=1.0/valA[j];
                csrColIndA[i]=i;
                csrRowIndA[i]=i;
            }
        }
    }
    csrRowIndA[num_rows]=nnz;

    myTopo=_A->getMyTopo();
    myArch=_A->getMyArch();

}
