#pragma once

#include <Eigen/Dense>

class Particle{

    private:
    //static int _count;

    public:

    Eigen::Vector3d pos;    //Position of charge
    Eigen::Vector3d com;    //COM position
    Eigen::Vector3d qDisp;  //Charge vector
    std::string name;       //name
    double r, b, q, rf;     //radius, length of charge vector, charge, minimum distance to wall

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

        this->pos = this->qDisp + this->com;
    };


    void chargeTrans(double step){
        Eigen::Vector3d v = Random::get_norm_vector();

        this->qDisp += step * v;
        if(this->qDisp.norm() > this->r){
            this->qDisp = this->qDisp.normalized() * this->r;
        }
        this->b = this->qDisp.norm();

        this->pos = this->qDisp + this->com;
    };


    void chargeTransRand(){
        //printf("Particle Move\n");
        this->qDisp = Random::get_norm_vector();
        //printf("raw: %lf\n", this->qDisp.norm());
        this->qDisp = this->qDisp * Random::get_random() * this->r;
        //printf("fixed: %lf\n", this->qDisp.norm());
        this->b = this->qDisp.norm();
        this->pos = this->qDisp + this->com;
        //printf("particle Move end\n");
    };


    void translate(std::vector<double> v){
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
