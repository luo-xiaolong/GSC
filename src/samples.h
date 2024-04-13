#pragma once

#include <iostream>
#include <map>
#include <cstring>
#include <string.h>
#include <vector>
using namespace std;


class Samples{
     std::map<std::string, uint32_t> whichIndMap;
     bool samples_init;
     string all_samples;
public:

    Samples(){
        samples_init=false;
        all_samples="";
		no_samples=0;
    }
    uint32_t no_samples;
    
    int loadSamples(vector<string>& v_samples);  
    uint32_t getWhich(std::string nm);

    uint32_t * setSamples(const std::string  & samples, string &str);
    void get_all_samples(string &str);

};



