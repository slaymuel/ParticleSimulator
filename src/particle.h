#pragma once

#include "random/random.h"

namespace Simulator{
class Particle{

    private:

    public:

    // Position of charge
    Eigen::Vector3d pos;    
    // center of mass position
    Eigen::Vector3d com;  
    // Local position of the charge 
    Eigen::Vector3d qDisp;  
    // Name of the particle
    std::string name;       
    //radius, length of charge vector, charge, minimum distance to wall
    double r, b, b_min, b_max, q, rf;     
    // The index of the particle in the particles object
    unsigned int index;

    // Translate the particle in a random direction
    // with a magnitude between 0 and step
    void translate(double step);
    // Translate the particle charge in a random
    // direction with maximum magnitude step
    void chargeTranslate(double step);
    // Moves the particle charge to a random position
    // which lies inside the bounds of the particle
    void chargeTransRand();
    // Rotate the particle
    void rotate(double step);
    // Translate the particle according to v
    template<typename T>
    void translate(T &v){
        this->com[0] +=  v[0];
        this->com[1] +=  v[1];
        this->com[2] +=  v[2];
        this->pos = this->qDisp + this->com;
    }

};

}