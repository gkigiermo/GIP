#ifndef _gip_timercuda_
#define _gip_timercuda_

#include <sys/time.h>
#include <sys/times.h>
#include "../../AbstractInterfaces/GIP_Timer.h"
#include "Headers/GIP_gpu.h"

#include "cuda.h"
#include "cuda_runtime.h"


class GIP_TimerCuda : public GIP_Timer{
    public:
        GIP_TimerCuda();
        void startChannel();
        void stopChannel();
        void stopAllChannels();

    protected:
	cudaEvent_t start_t;
	cudaEvent_t stop_t;
	
};

#endif
