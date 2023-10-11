#include "GIP_Timer.h"


GIP_Timer::GIP_Timer()
{
    tm_l=0;
    sum_l=0;
    tm_g=0;
    sum_g=0;
    cont=0;
}

float GIP_Timer::getTimeL()
{
    return tm_l;
}
float GIP_Timer::getSumL()
{
    return sum_l;
}
float GIP_Timer::getTimeG()
{
    return tm_g;
}
float GIP_Timer::getSumG()
{
    return sum_g;
}
float GIP_Timer::getTimeAvg()
{
    if(cont!=0)
    {
        return sum_g/cont;
    }
    else return 0;
}
