#pragma once

#include "random.h"

using Vector3 = Eigen::Block<Eigen::MatrixXd, 1, 3>;

class Particle{

    private:
    //static int _count;

    public:

    Eigen::Vector3d pos;    //Position of charge
    Eigen::Vector3d com;    //COM position
    Eigen::Vector3d qDisp;  //Charge vector
    std::string name;       //name
    double r, b, b_min, b_max, q, rf;     //radius, length of charge vector, charge, minimum distance to wall

    unsigned int index;

    /*Particle& operator=(const Particle& p){
        this->pos = p.pos;
        this->r = p.r;
        this->rf = p.rf;
        this->q = p.q;
        this->b = p.b;
        this->index = p.index;

        return *this;
    }*/


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
        /*Eigen::Vector3d v = Random::get_norm_vector();
        this->qDisp += step * v * Random::get_random();

        this->b = this->qDisp.norm();

        this->pos = this->qDisp + this->com;*/

        Eigen::Vector3d v = Random::get_vector();
        this->pos[0] += step * v[0];
        this->pos[1] += step * v[1];
        this->pos[2] += step * v[2];
    };


    void chargeTransRand(){
        
        //this->b = this->b_min + (this->b_max - this->b_min) * Random::get_random();
        //this->qDisp = Random::get_random_vector(1.0);
        //this->qDisp = this->qDisp.normalized() * this->b;

        //this->qDisp = Random::get_norm_vector();
        //this->qDisp *= (this->b_min + (this->b_max - this->b_min) * Random::get_random());
        
        this->qDisp = Random::get_random_vector(this->b_max);
        this->b = this->qDisp.norm();

        //this->qDisp = Random::get_random_vector(0.5);
        //this->qDisp *= 2.0 * this->b_max;
        //this->b = this->qDisp.norm();

        this->pos = this->com + this->qDisp;
    };


    void translate(std::vector<double> &v){
        this->com[0] +=  v[0];
        this->com[1] +=  v[1];
        this->com[2] +=  v[2];
        this->pos = this->qDisp + this->com;
    }

    void translate(Eigen::Vector3d &v){
        this->com +=  v;
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
