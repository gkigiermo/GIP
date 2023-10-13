#include "cuda.h"
#include "cuda_runtime.h"

extern "C++" void GIP_cudaMallocDouble(double* &,int );

extern "C++" void GIP_cudaMallocInt(int* &,int );

extern "C++" void GIP_cudaFree(void* );

extern "C++" void GIP_cudaMemcpyDToGpu(double* ,double* ,int );

extern "C++" void GIP_cudaMemcpyIToGpu(int* ,int* ,int );

extern "C++" void GIP_cudaMemcpyIToCpu(int* , int* ,int );

extern "C++" void GIP_cudaMemcpyDToCpu(double* , double* ,int );

extern "C++" void GIP_cudaGetLastError();

extern "C++" void GIP_cudaDcopy(int, double*, double *);

extern "C++" void GIP_cudaDCsrSpMV(int , int *, int* , double* , double* , double* ,int ,int);

extern "C++" void GIP_cudaDCsrSpMV(int , int *, int* , double* , double* , double* ,double,double,int ,int);

extern "C++" void GIP_cudaDaxpy(int ,double* ,double * ,double ,int ,int );

extern "C++" void GIP_cudaDaxpy(int ,double* ,double * ,double ,double ,int ,int );

extern "C++" double GIP_cudaDdot(int , double* ,double* ,double* ,double* , int ,int );

extern "C++" void  GIP_cudaEventCreate(cudaEvent_t&);

extern "C++" void  GIP_cudaEventRecord(cudaEvent_t&);

extern "C++" void  GIP_cudaEventSynchronize(cudaEvent_t& );

extern "C++" float GIP_cudaEventElapsedTime(cudaEvent_t& ,cudaEvent_t&);

extern "C++" void  GIP_cudaEventDestroy(cudaEvent_t& );

extern "C++" void  GIP_cudaSetDevice(int );

extern "C++" double GIP_gpuNormi(double* ,double* ,double* ,int , int ,int );

extern "C++" double GIP_gpuGetCFL(double* ,double* ,double* ,double* ,double* ,double* ,int ,double ,double , int ,int );

