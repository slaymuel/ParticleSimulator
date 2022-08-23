#include "particle.h"

namespace Simulator{

void Particle::translate(double step){
    Eigen::Vector3d v = Random::get_vector();
    this->com[0] += step * v[0];
    this->com[1] += step * v[1];
    this->com[2] += step * v[2];

    this->pos = this->com + this->qDisp;
}

void Particle::chargeTranslate(double step){
    Eigen::Vector3d v = Random::get_vector();
    this->pos[0] += step * v[0];
    this->pos[1] += step * v[1];
    this->pos[2] += step * v[2];
}


void Particle::chargeTransRand(){
    this->qDisp = Random::get_random_vector(this->b_max);
    this->b = this->qDisp.norm();

    this->pos = this->com + this->qDisp;
}

void Particle::rotate(double step){
    Eigen::Vector3d v = Random::get_vector();
    //v *= step;
    this->qDisp += v * step;
    this->qDisp = this->qDisp.normalized() * this->b;
    this->pos = this->com + this->qDisp;
}


}