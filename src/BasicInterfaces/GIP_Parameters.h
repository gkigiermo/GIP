#ifndef  _gip_parameters_
#define  _gip_parameters_
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
using namespace std;

class GIP_Parameters
{
    friend ostream &operator<<(ostream &, const GIP_Parameters &);
    
  public:
    GIP_Parameters();
    GIP_Parameters(const GIP_Parameters&);
    const GIP_Parameters& operator=(const GIP_Parameters&);
    ~GIP_Parameters();

    inline const bool       getBool(const string& key) const;
    inline const int        getInt(const string& key) const;
    inline const double&     getDouble(const string& key) const;    
    inline const string&     getString(const string& key) const; 
    inline const char*       getChar(const string& key) const; 

    void insertBool(const string&key, const bool value);
    void insertInt(const string& key, const int value);
    void insertDouble(const string& key, const double& value);
    void insertString(const string& key, const string& value);

    void setBool(const string& key, const bool value);
    void setInt(const string& key, const int value);
    void setDouble(const string& key, const double& value);
    void setString(const string& key, const string& value);


    protected:

    map<string, bool>      bools;
    map<string, int>       integers;
    map<string, double>    doubles;
    map<string,string>     strings;

    string name;
};

//inlines
inline const bool GIP_Parameters::getBool(const string& key) const{
    map<string,bool>::const_iterator i=bools.find(key);
    return i->second;
}
inline const int GIP_Parameters::getInt(const string& key) const{
    map<string,int>::const_iterator i=integers.find(key);
    return i->second;
}
inline const double&  GIP_Parameters::getDouble(const string& key) const{
    map<string,double>::const_iterator i=doubles.find(key);
    return i->second;
}
inline const string&  GIP_Parameters::getString(const string& key) const{
    map<string,string>::const_iterator i=strings.find(key);
    return i->second;
}
inline const char*  GIP_Parameters::getChar(const string& key) const{
    map<string,string>::const_iterator i=strings.find(key);
    return i->second.data();
}

#endif
