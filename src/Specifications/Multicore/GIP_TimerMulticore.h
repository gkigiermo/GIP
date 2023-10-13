#ifndef _gip_timermulticore_
#define _gip_timermulticore_
#include <stdlib.h>
#include "../../AbstractInterfaces/GIP_Timer.h"


class GIP_TimerMulticore : public GIP_Timer{
    public:
        GIP_TimerMulticore();
        void startChannel();
        void stopChannel();
        void stopAllChannels();

    protected:

        double getTime();
        double end;
        double begin;
};

#endif
