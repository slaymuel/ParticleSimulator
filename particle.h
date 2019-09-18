#pragma once

#include <Eigen/Dense>

class Particle{

    private:
    //static int _count;

    public:

    Eigen::Vector3d pos;
    Eigen::Vector3d com;
    std::string name;
    double r, b, q, rf;

    unsigned int index;

    Particle& operator=(const Particle& p){
        this->pos = p.pos;
        this->r = p.r;
        this->rf = p.rf;
        this->q = p.q;
        this->b = p.b;
        this->index = p.index;

        return *this;
    }

    void translate(double step){
        Eigen::Vector3d v = Random::get_vector();
        this->pos[0] += step * (v[0] * 2.0 - 1);
        this->pos[1] += step * (v[1] * 2.0 - 1);
        this->pos[2] += step * (v[2] * 2.0 - 1);
    };

    void rotate(double step){
        printf("Rotating particle %i\n", this->index);
    };
};
