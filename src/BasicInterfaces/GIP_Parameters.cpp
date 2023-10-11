#include "GIP_Parameters.h"

GIP_Parameters::GIP_Parameters(){
    name="";
}

GIP_Parameters::GIP_Parameters(const GIP_Parameters& o) : bools(o.bools), integers(o.integers), doubles(o.doubles), strings(o.strings),name(o.name){}
    
const GIP_Parameters& GIP_Parameters::operator=(const GIP_Parameters& o){
    bools=o.bools;
    integers=o.integers;
    doubles=o.doubles; 
    strings=o.strings;
    name=o.name;    
    return *this;
}

GIP_Parameters::~GIP_Parameters(){ 
    map<string, bool>().swap(bools);
    map<string, int>().swap(integers);
    map<string, double>().swap(doubles);
    map<string,string>().swap(strings);
}


void GIP_Parameters::insertBool(const string& key,const bool value){
    bools[key]=value;
}
void GIP_Parameters::insertInt(const string& key,const int value){
    integers[key]=value;
}
void GIP_Parameters::insertDouble(const string& key,const double& value){
    doubles[key]=value;    
}
void GIP_Parameters::insertString(const string& key,const string& value){
    strings[key]=value;    
}
void GIP_Parameters::setBool(const string& key,const bool value){
    bools[key]=value;
}
void GIP_Parameters::setInt(const string& key,const int value){
    integers[key]=value;
}
void GIP_Parameters::setDouble(const string& key,const double& value){
    doubles[key]=value;    
}
void GIP_Parameters::setString(const string& key,const string& value){
    strings[key]=value;    
}

ostream &operator<<(ostream & qout, const GIP_Parameters& p){
    qout<<"   START PARAMETERS "<<p.name<<endl<<endl;
    qout<<"    -Booleans:"<<endl;
    for(map<string,bool>::const_iterator i=p.bools.begin(); i!=p.bools.end(); ++i)
        qout<<"     "<<i->first<<": "<<i->second<<endl;
    qout<<endl<<"    -Integers:"<<endl;
    for(map<string,int>::const_iterator i=p.integers.begin(); i!=p.integers.end(); ++i)
        qout<<"     "<<i->first<<": "<<i->second<<endl;                                          //aqui no podrÃ© ficar els glovals ja que no els tinc
    qout<<endl<<"    -Doubles:"<<endl;
    for(map<string,double>::const_iterator i=p.doubles.begin(); i!=p.doubles.end(); ++i)
        qout<<"     "<<i->first<<": "<<i->second<<endl; //aqui no podrÃ© ficar els glovals ja que no els tinc
    qout<<endl<<"    -Strings:"<<endl;
    for(map<string,string>::const_iterator i=p.strings.begin(); i!=p.strings.end(); ++i)
        qout<<"     "<<i->first<<": "<<i->second<<endl; 
    qout<<endl<<"   END PARAMETERS"<<endl;
    return qout;
}
