#ifndef _gip_timer_
#define _gip_timer_

#include <mpi.h>

class GIP_Timer{
    public:

        GIP_Timer();
        virtual void startChannel(){};
        virtual void stopChannel(){};
        virtual void stopAllChannels(){};
        float getTimeL();
        float getSumL();
        float getTimeG();
        float getSumG();
        float getTimeAvg();

    protected:

        float tm_l;
        float sum_l;
        float tm_g;
        float sum_g;
        int cont;
};

#endif
