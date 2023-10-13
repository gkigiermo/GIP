#include "GIP_MatrixCuda.h"

GIP_MatrixCuda::GIP_MatrixCuda(char* name,GIP_Arch* _architecture): GIP_Matrix(name,_architecture) {
	GIP_cudaMallocDouble(dcsrValA,nnz);
	GIP_cudaMallocInt(dcsrColIndA,nnz);
	GIP_cudaMallocInt(dcsrRowIndA,num_rows+1);
	GIP_cudaMemcpyDToGpu(csrValA,dcsrValA,nnz);
	GIP_cudaMemcpyIToGpu(csrColIndA,dcsrColIndA,nnz);
	GIP_cudaMemcpyIToGpu(csrRowIndA,dcsrRowIndA,num_rows+1);
	threads=128;
	blocks=(num_rows+(threads-1))/threads;
}

void GIP_MatrixCuda::postConstruct(string matrix_name,GIP_Topo* _topo,GIP_Arch* _architecture)
{
	int rank,nz;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nz);
	FILE *fp;

	fp= fopen(matrix_name.c_str(),"rb");

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

	GIP_cudaMallocDouble(dcsrValA,nnz);
	GIP_cudaMallocInt(dcsrColIndA,nnz);
	GIP_cudaMallocInt(dcsrRowIndA,num_rows+1);
	GIP_cudaMemcpyDToGpu(csrValA,dcsrValA,nnz);
	GIP_cudaMemcpyIToGpu(csrColIndA,dcsrColIndA,nnz);
	GIP_cudaMemcpyIToGpu(csrRowIndA,dcsrRowIndA,num_rows+1);
	threads=128;
	blocks=(num_rows+(threads-1))/threads;

	myTopo=_topo;
	myArch=_architecture;


}

GIP_MatrixCuda::~GIP_MatrixCuda()
{
    if(dcsrValA!=NULL)
        GIP_cudaFree(dcsrValA);
    if(dcsrColIndA!=NULL)
        GIP_cudaFree(dcsrColIndA);
    if(dcsrRowIndA!=NULL)
        GIP_cudaFree(dcsrRowIndA);
  
}

void GIP_MatrixCuda::printRows(int n)
{
    cout<<"Matrix Info Clase CUDA :"<<endl;
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


//Ax=b
void GIP_MatrixCuda::spmv(GIP_Vector* x, GIP_Vector* b)
{
    GIP_cudaDCsrSpMV(num_rows,dcsrRowIndA,dcsrColIndA,dcsrValA,x->getDevicePtr(),b->getDevicePtr(),blocks,threads);
}

void GIP_MatrixCuda::spmv(GIP_Vector* x, GIP_Vector* b,double alpha,double beta)
{
    GIP_cudaDCsrSpMV(num_rows,dcsrRowIndA,dcsrColIndA,dcsrValA,x->getDevicePtr(),b->getDevicePtr(),alpha,beta,blocks,threads);
}

void GIP_MatrixCuda::spmv(GIP_Vector* x, GIP_Matrix* Ab)
{
    GIP_cudaDCsrSpMV(num_rows,dcsrRowIndA,dcsrColIndA,dcsrValA,x->getDevicePtr(),Ab->getCsrValADevice(),blocks,threads);
}
void GIP_MatrixCuda::spmv(GIP_Vector* x, GIP_Matrix* Ab,double alpha,double beta)
{
    GIP_cudaDCsrSpMV(num_rows,dcsrRowIndA,dcsrColIndA,dcsrValA,x->getDevicePtr(),Ab->getCsrValADevice(),alpha,beta,blocks,threads);
}

void GIP_MatrixCuda::operator=(GIP_Matrix& matrix)
{
    num_rows=matrix.getNumRows();
    num_cols=matrix.getNumCols();
    nnz=     matrix.getNumNnz();

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

    GIP_cudaMallocDouble(dcsrValA,nnz);
    GIP_cudaMallocInt(dcsrColIndA,nnz);
    GIP_cudaMallocInt(dcsrRowIndA,num_rows+1);
    GIP_cudaMemcpyDToGpu(csrValA,dcsrValA,nnz);
    GIP_cudaMemcpyIToGpu(csrColIndA,dcsrColIndA,nnz);
    GIP_cudaMemcpyIToGpu(csrRowIndA,dcsrRowIndA,num_rows+1);
    threads=128;
    blocks=(num_rows+(threads-1))/threads;


    myTopo=matrix.getMyTopo();
    myArch=matrix.getMyArch();

}

void GIP_MatrixCuda::setUpDiagonal(GIP_Matrix* _A)
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

    GIP_cudaMallocDouble(dcsrValA,nnz);
    GIP_cudaMallocInt(dcsrColIndA,nnz);
    GIP_cudaMallocInt(dcsrRowIndA,num_rows+1);
    GIP_cudaMemcpyDToGpu(csrValA,dcsrValA,nnz);
    GIP_cudaMemcpyIToGpu(csrColIndA,dcsrColIndA,nnz);
    GIP_cudaMemcpyIToGpu(csrRowIndA,dcsrRowIndA,num_rows+1);
    threads=128;
    blocks=(num_rows+(threads-1))/threads;

    myTopo=_A->getMyTopo();
    myArch=_A->getMyArch();
}
