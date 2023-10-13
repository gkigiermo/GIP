#ifndef _gip_arch_
#define _gip_arch_

class GIP_Arch{

    public:
        GIP_Arch(){};
        GIP_Arch(int, int, int);
        ~GIP_Arch(){};
        virtual void postConstruct(int,int,int){};
        virtual void setUp(){};

    protected:
        int num_cpus;
        int num_cores;
        int num_devs;
};

#endif
