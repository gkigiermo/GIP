#include<iostream>
#include<cuda.h>
#include<cuda_runtime.h>
#include<cuda_runtime_api.h>
#include<cublas.h>

using namespace std;
/*
extern "C++" void GIP_cusparseDcsrmv(double *d_val,int *d_indxCol,int *d_indxRow,double *d_x,double *d_y,int rows,int cols,int nnz,cusparseMatDescr_t *descr,cusparseHandle_t *handle)
{
        double alpha=1.0;
        double beta=0.0;
       cusparseDcsrmv(*handle, CUSPARSE_OPERATION_NON_TRANSPOSE, rows, cols,nnz, &alpha,*descr, d_val, d_indxRow, d_indxCol, d_x, &beta, d_y);
}
*/

// Para reservar memoria Double e Int
extern "C++" void GIP_cudaMallocDouble(double* &vector,int size)
{        
	cudaMalloc((void**)&vector,size*sizeof(double));
}

extern "C++" void GIP_cudaMallocInt(int* &vector,int size)
{        
	cudaMalloc((void**)&vector,size*sizeof(int));
}

// Para copiar a CPU->GPU Double e Int
extern "C++" void GIP_cudaMemcpyDToGpu(double* h_vect,double* d_vect,int size )
{
		cudaMemcpy(d_vect,h_vect,size*sizeof(double),cudaMemcpyHostToDevice);
		
}

extern "C++" void GIP_cudaMemcpyIToGpu(int* h_vect,int* d_vect,int size )
{
		cudaMemcpy(d_vect,h_vect,size*sizeof(int),cudaMemcpyHostToDevice);
		
}
// Para copiar a GPU->CPU Double e Int
extern "C++" void GIP_cudaMemcpyIToCpu(int* h_vect, int* d_vect,int size )
{
		cudaMemcpy(h_vect,d_vect,size*sizeof(int),cudaMemcpyDeviceToHost);
}

extern "C++" void GIP_cudaMemcpyDToCpu(double* h_vect, double* d_vect,int size )
{
                cudaMemcpy(h_vect,d_vect,size*sizeof(double),cudaMemcpyDeviceToHost);
}

// Para liberar memoria
extern "C++" void GIP_cudaFree(void* vector)
{
	cudaFree(vector);
}

extern "C++" void GIP_cudaGetLastError(){
     cudaError_t error;
     error=cudaGetLastError();
     if(error!= cudaSuccess)
     {
       cout<<" ERROR INSIDE A CUDA FUNCTION: "<<error<<" "<<cudaGetErrorString(error)<<endl;
       exit(0);
     }
}

//vec1=vec2
extern "C++" void GIP_cudaDcopy(int size, double* vec1, double* vec2)
{
	cudaMemcpy(vec1,vec2,size*sizeof(double),cudaMemcpyDeviceToDevice);
}

/* In house implementations */
// y=Ax
__global__ void cudaDcsrspmv(int num_rows, int *rowIndA, int* colIndA,double* valA,double* x, double* y)
{
    int row =blockDim.x*blockIdx.x +threadIdx.x;
    if( row<num_rows)
    {
        double sum =0;
        int row_start= rowIndA[row];
        int row_end= rowIndA[row+1];
        for(int j=row_start;j<row_end;j++)
            sum+=valA[j]*x[colIndA[j]];
        y[row]=sum;
    }
}

extern "C++" void GIP_cudaDCsrSpMV(int num_rows, int *rowIndA, int* colIndA,double* valA,double* x, double* y,int blocks,int threads)
{
    dim3 dimGrid(blocks,1,1);
    dim3 dimBlock(threads,1,1);
    cudaFuncSetCacheConfig( cudaDcsrspmv, cudaFuncCachePreferL1 ); // para asignar 48KB a cache en el kernel csr
    cudaDcsrspmv<<<dimGrid,dimBlock>>>(num_rows,rowIndA,colIndA,valA,x,y);

}

// y=alpha*Ax+beta*y
__global__ void cudaDcsrspmvab(int num_rows, int *rowIndA, int* colIndA,double* valA,double* x, double* y, double alpha,double beta)
{
    int row =blockDim.x*blockIdx.x +threadIdx.x;
    if( row<num_rows)
    {
        double sum =0;
        int row_start= rowIndA[row];
        int row_end= rowIndA[row+1];
        for(int j=row_start;j<row_end;j++)
            sum+=valA[j]*x[colIndA[j]];
        y[row]=alpha*sum+beta*y[row];
    }
}

extern "C++" void GIP_cudaDCsrSpMV(int num_rows, int *rowIndA, int* colIndA,double* valA,double* x, double* y,double alpha,double beta,int blocks,int threads)
{
    dim3 dimGrid(blocks,1,1);
    dim3 dimBlock(threads,1,1);
    cudaFuncSetCacheConfig( cudaDcsrspmvab, cudaFuncCachePreferL1 ); // para asignar 48KB a cache en el kernel csr
    cudaDcsrspmvab<<<dimGrid,dimBlock>>>(num_rows,rowIndA,colIndA,valA,x,y,alpha,beta);

}


// y=ax+y
__global__ void cudaDaxpy(int n, double* vec1, double *vec2, double alpha)
{
    int row =blockDim.x*blockIdx.x +threadIdx.x;
    if(row<n)
    {
        vec1[row]+=alpha*vec2[row];
    }
}

extern "C++" void GIP_cudaDaxpy(int n, double* vec1, double *vec2,double alpha,int blocks, int threads)
{

    dim3 dimGrid(blocks,1,1);
    dim3 dimBlock(threads,1,1);
    cudaDaxpy<<<dimGrid,dimBlock>>>(n,vec1,vec2,alpha);
}


//y=ax+by
__global__ void cudaDaxpby(int n, double* vec1, double *vec2, double alpha,double beta)
{
    int row =blockDim.x*blockIdx.x +threadIdx.x;
    if(row<n)
    {
        vec1[row]=alpha*vec2[row]+beta*vec1[row];
    }
}

extern "C++" void GIP_cudaDaxpy(int n, double* vec1, double *vec2,double alpha,double beta,int blocks, int threads)
{

    dim3 dimGrid(blocks,1,1);
    dim3 dimBlock(threads,1,1);
   cudaDaxpby<<<dimGrid,dimBlock>>>(n,vec1,vec2,alpha,beta);
}

__global__ void cudaDdot(double *g_idata1, double *g_idata2, double *g_odata, unsigned int n)
{
    extern __shared__ double sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    double mySum = (i < n) ? g_idata1[i]*g_idata2[i] : 0;


    if (i + blockDim.x < n)
        mySum += g_idata1[i+blockDim.x]*g_idata2[i+blockDim.x];
    
    sdata[tid] = mySum;
    __syncthreads();

    for (unsigned int s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
            sdata[tid] = mySum = mySum + sdata[tid + s];

        __syncthreads();
    }

    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


extern "C++" double GIP_cudaDdot(int n, double* vec1,double* vec2,double* h_temp,double* d_temp, int blocks,int threads)
{

    dim3 dimGrid(blocks,1,1);
    dim3 dimBlock(threads,1,1);
    cudaDdot<<<dimGrid,dimBlock,threads*sizeof(double)>>>(vec1,vec2,d_temp,n);


    cudaMemcpy(h_temp, d_temp, blocks * sizeof(double), cudaMemcpyDeviceToHost);

    double sum=0;
    for(int i=0;i<blocks;i++)
    {
        sum+=h_temp[i];
    }
    return sum;
}


extern "C++" void GIP_cudaEventCreate(cudaEvent_t& evet)
{
	cudaEventCreate(&evet);
}

extern "C++" void GIP_cudaEventRecord(cudaEvent_t& evet)
{
	cudaEventRecord(evet,0);
}

extern "C++" void GIP_cudaEventRecord(cudaEvent_t& evet, cudaStream_t& streamt)
{
	cudaEventRecord(evet,streamt);
}

extern "C++" void GIP_cudaEventSynchronize(cudaEvent_t& evet)
{
	cudaEventSynchronize(evet);
}

extern "C++" float GIP_cudaEventElapsedTime(cudaEvent_t& start,cudaEvent_t &stop)
{
	float tpo;
	cudaEventElapsedTime(&tpo,start,stop);
	return tpo;
}

extern "C++" void GIP_cudaEventDestroy(cudaEvent_t& evet)
{
	cudaEventDestroy(evet);
}

extern "C++" void GIP_cudaSetDevice(int numDev)
{
	cudaSetDevice(numDev);
}


__global__ void getCFL(double *u, double* v,double *w,double *dxs,double *g_odata, unsigned int n,double gamma, double rho)
{
    // double *sdata = SharedMemory<double>();
    extern __shared__ double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    double dt=1e10;

    double convCoef = 0.25;

    double up=(i<n) ? fabs(u[i]) :0;
    double vp=(i<n) ? fabs(v[i]) :0;
    double wp=(i<n) ? fabs(w[i]) :0;
    double dx=(i<n) ? dxs[i] :0;

    double normiVel= (up>vp) ? up : vp;
    normiVel= (wp>normiVel) ? wp : normiVel;
    if(!(dx/normiVel < 1e-17))
    {
        if(dt>convCoef*dx/(normiVel))
            dt= convCoef*dx/(normiVel);
    }
    if(!(dx*dx*rho/gamma  < 1e-17))
    {
        if(dt>0.2*dx*dx*rho/gamma)
            dt=0.2*dx*dx*rho/gamma;
    }

    if (i + blockDim.x < n)
    {

        up=fabs(u[i+blockDim.x]) ;
        vp=fabs(v[i+blockDim.x]) ;
        wp=fabs(w[i+blockDim.x]) ;
        dx=dxs[i+blockDim.x] ;

        normiVel= (up>vp) ? up : vp;
        normiVel= (wp>normiVel) ? wp : normiVel;
        if(!(dx/normiVel < 1e-17))
        {
            if(dt>convCoef*dx/(normiVel))
                dt= convCoef*dx/(normiVel);
        }
        if(!(dx*dx*rho/gamma  < 1e-17))
        {
            if(dt>0.2*dx*dx*rho/gamma)
                dt=0.2*dx*dx*rho/gamma;
        }

    } 
    sdata[tid] = dt;
    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
           sdata[tid] = dt = (sdata[tid + s]<dt) ? sdata[tid + s] : dt;
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


extern "C++" double GIP_gpuGetCFL(double* d_u,double* d_v,double* d_w,double* d_dxs,double* h_odata,double* d_odata,int n,double gamma,double rho, int threads,int blocks)
{

            dim3 dimGrid(blocks,1,1);
            dim3 dimBlock(threads,1,1);// (1024,1024,64)
            getCFL<<<dimGrid,dimBlock,threads*sizeof(double)>>>(d_u,d_v,d_w,d_dxs,d_odata,n,gamma,rho);

           
            cudaMemcpy(h_odata, d_odata, blocks * sizeof(double), cudaMemcpyDeviceToHost);    

            double min=1e10;
            for(int i=0;i<blocks;i++)
            {
                     min=(min>h_odata[i]) ? h_odata[i] : min;
            }

            
            min=(0.8*min>1e-1)? 1e-1 :0.8*min;
            return min;
}


__global__ void normi(double *g_idata, double *g_odata, unsigned int n)
{
   // double *sdata = SharedMemory<double>();
    extern __shared__ double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    double myMax = (i < n) ? fabs(g_idata[i]) : 0;

    if (i + blockDim.x < n)
    {
        double next =fabs( g_idata[i+blockDim.x]);
        myMax=(next>myMax) ? next : myMax;
    }
    // mySum += g_idata[i+blockDim.x];

    sdata[tid] = myMax;
    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {

           sdata[tid] = myMax = (sdata[tid + s]>myMax) ? sdata[tid + s] : myMax;

           // sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


 extern "C++" double GIP_gpuNormi(double* d_val,double* h_odata,double* d_odata,int n, int threads,int blocks)
{

            dim3 dimGrid(blocks,1,1);
            dim3 dimBlock(threads,1,1);// (1024,1024,64)
            normi<<<dimGrid,dimBlock,threads*sizeof(double)>>>(d_val,d_odata,n);

           
            cudaMemcpy(h_odata, d_odata, blocks * sizeof(double), cudaMemcpyDeviceToHost);    

            double max=0;
            for(int i=0;i<blocks;i++)
            {
                     max=(max<h_odata[i]) ? h_odata[i] : max;
            }
            return max;
}
