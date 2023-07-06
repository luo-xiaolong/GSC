
#include "defs.h"
#include "samples.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>
// *****************************************************************************************************************
int Samples::loadSamples(vector<string>& v_samples)  
{
    uint32_t c = 0;
    for (size_t i = 0; i < v_samples.size(); i++){
        if (!whichIndMap.insert(std::pair<std::string, uint32_t>(v_samples[i], c++)).second) {
                std::cout << "Error! Two individuals with the same name!\n";
                return 2;
            }
        all_samples  += v_samples[i] + "\t";
    }            
	all_samples  =  all_samples.substr(0, all_samples.length() -1);    
   
    no_samples = v_samples.size();
    return 0;
}
// *****************************************************************************************************************
void Samples::get_all_samples(string &str){
    str+=all_samples;
}
// *****************************************************************************************************************
uint32_t Samples::getWhich(std::string nm)
{
    auto it = whichIndMap.find(nm);
    if(it != whichIndMap.end())
             return it->second;
    else
    {
        std::cout << "There is no sample " << nm << " in the set\n";
        exit(1);
        
    }
}
// *****************************************************************************************************************
uint32_t * Samples::setSamples(const std::string & samples, string &str)
{
    uint32_t * smplIDs = nullptr;
    long size = 0;
    
    no_samples = 0;
    if(samples[0] == '@')
    {
        std::ifstream in_samples(samples.substr(1));
        if(!in_samples.is_open())
        {
            std::cout << "Error. Cannot open " << samples.substr(1)<< " file with samples.\n";
            exit(1);
        }
        
        uint32_t i = 0;
        std::string item;
        
        while (in_samples >> item) {
            size++;
            
        }
        cout<<"size:"<<size<<endl;
        in_samples.clear();
        in_samples.seekg(0);
        
        smplIDs = new uint32_t[size];
        while (in_samples >> item) {
            str+=(item+"\t");
            smplIDs[i++] = getWhich(item);
            no_samples++;
        }
    }
    else
    {
        size = std::count(samples.begin(), samples.end(), ',') + 1;
        
        smplIDs = new uint32_t[size];
        char delim = ',';
        uint32_t i = 0;
        std::stringstream ss(samples);
        std::string item;
        while (getline(ss, item, delim)) {
            str+=(item+"\t");
            smplIDs[i++] = getWhich(item);
            no_samples++;
        }
    }
	str = str.substr(0, str.length() -1);
    
    return smplIDs;
}
