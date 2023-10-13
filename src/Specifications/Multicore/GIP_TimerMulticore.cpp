#include "GIP_TimerMulticore.h"

#include <omp.h>

GIP_TimerMulticore::GIP_TimerMulticore() : GIP_Timer() {}

void GIP_TimerMulticore::startChannel()
{
    begin=getTime();
}
void GIP_TimerMulticore::stopChannel()
{
    end=getTime();
    tm_l=(float)(end-begin);      
    sum_l=sum_l+tm_l;
    cont=cont+1;
}
void GIP_TimerMulticore::stopAllChannels()
{
    end=getTime();

    tm_l=(float)(end-begin); 
    sum_l=sum_l+tm_l;
    MPI_Allreduce(&tm_l,&tm_g,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);    
    sum_g=sum_g+tm_g;
    cont=cont+1;
}

double GIP_TimerMulticore::getTime()
{
    return omp_get_wtime();
}

