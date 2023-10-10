#include "GIP_TimerCuda.h"

GIP_TimerCuda::GIP_TimerCuda() : GIP_Timer() {
	GIP_cudaEventCreate(start_t);
	GIP_cudaEventCreate(stop_t);
}

void GIP_TimerCuda::startChannel()
{
	GIP_cudaEventRecord(start_t);
}
void GIP_TimerCuda::stopChannel()
{
	GIP_cudaEventRecord(stop_t);
	GIP_cudaEventSynchronize(stop_t);
	
	tm_l=GIP_cudaEventElapsedTime(start_t,stop_t)/1000.0;
	sum_l=sum_l+tm_l;
	cont=cont+1;

}
void GIP_TimerCuda::stopAllChannels()
{
	GIP_cudaEventRecord(stop_t);
	GIP_cudaEventSynchronize(stop_t);
	
	tm_l=GIP_cudaEventElapsedTime(start_t,stop_t)/1000.0;
	sum_l=sum_l+tm_l;
	MPI_Allreduce(&tm_l,&tm_g,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);	
	sum_g=sum_g+tm_g;
	cont=cont+1;
}

