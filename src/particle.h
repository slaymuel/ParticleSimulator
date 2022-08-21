#pragma once

#include "random/random.h"

namespace Simulator{

using Vector3 = Eigen::Block<Eigen::MatrixXd, 1, 3>;

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

    unsigned int index;


    void translate(double step){
        Eigen::Vector3d v = Random::get_vector();
        this->com[0] += step * v[0];
        this->com[1] += step * v[1];
        this->com[2] += step * v[2];

        this->pos = this->com + this->qDisp;
    };

    void chargeTranslate(double step){
        Eigen::Vector3d v = Random::get_vector();
        this->pos[0] += step * v[0];
        this->pos[1] += step * v[1];
        this->pos[2] += step * v[2];
    };

    void chargeTrans(double step){
        Eigen::Vector3d v = Random::get_vector();
        this->pos[0] += step * v[0];
        this->pos[1] += step * v[1];
        this->pos[2] += step * v[2];
    };


    void chargeTransRand(){
        this->qDisp = Random::get_random_vector(this->b_max);
        this->b = this->qDisp.norm();

        this->pos = this->com + this->qDisp;
    };

    template<typename T>
    void translate(T &v){
        this->com[0] +=  v[0];
        this->com[1] +=  v[1];
        this->com[2] +=  v[2];
        this->pos = this->qDisp + this->com;
    }

    void rotate(double step){
        Eigen::Vector3d v = Random::get_vector();
        //v *= step;
        this->qDisp += v * step;
        this->qDisp = this->qDisp.normalized() * this->b;
        this->pos = this->com + this->qDisp;
    };
};

}