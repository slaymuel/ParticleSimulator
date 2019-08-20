#pragma once

#include <Eigen/Dense>

class Particle{

    private:
    //static int _count;

    public:

    Eigen::VectorXd pos;
    Eigen::VectorXd* chargePos;
    double r;
    double q;
    double b;
    int index;

    Particle& operator=(const Particle& p){
        this->pos = p.pos;
        this->r = p.r;
        this->q = p.q;
        this->b = p.b;
        this->index = p.index;

        return *this;
    }

    void translate(){
        this->pos += Random::get_random_vector();
        printf("Translating particle %i\n", this->index);
    };

    void rotate(){
        printf("Rotating particle %i\n", this->index);
    };

    double distance(std::shared_ptr<Particle> other){
        return std::sqrt( (this->pos[0] - other->pos[0]) * (this->pos[0] - other->pos[0])
                         + (this->pos[1] - other->pos[1]) * (this->pos[1] - other->pos[1]) 
                         + (this->pos[2] - other->pos[2]) * (this->pos[2] - other->pos[2]) );
    }
};
