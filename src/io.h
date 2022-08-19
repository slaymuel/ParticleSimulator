#pragma once

#include <iomanip>
#include <fstream>

#include "particles.h"

namespace Simulator{
    
namespace IO{

    void to_gro(std::string fileName, const Particles& p, const std::vector<double> d);

    void to_cpt(std::string fileName, const Particles& p, const std::vector<double> d);

    template<typename T>
    void to_xyz(std::string fileName, const T& p, const std::vector<double> d){
        // Write charge positions
        std::ofstream f (fileName + ".xyz");
        if (f.is_open())
        {
            f << p.size() << "\n\n";
            for(unsigned int i = 0; i < p.size(); i++)
                f << std::fixed << std::setprecision(3) << p[i]->name << " " <<  p[i]->pos[0] << " " << 
                                                           p[i]->pos[1] << " " << p[i]->pos[2] << "\n";
            f << "10 10 10" << "\n";
            f.close();
        }
        else std::cout << "Unable to open file\n";

        // Write center of mass positions
        std::ofstream fq (fileName + "_com.xyz");
        if (fq.is_open())
        {
            fq << p.size() << "\n\n";
            for(unsigned int i = 0; i < p.size(); i++)
                fq << std::fixed << std::setprecision(3) << p[i]->name << " " <<  p[i]->com[0] << " " << 
                                                            p[i]->com[1] << " " << p[i]->com[2] << "\n";
            fq << d[0] << " " << d[1] << " " << d[2] << "\n";
            fq.close();
        }
        else std::cout << "Unable to open file\n";
    }
};

}