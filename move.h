#pragma once
#include "particle.h"
#include <math.h>

using CallBack = std::function<void(std::vector< std::shared_ptr<Particle> >)>;

class Move{
    public:

    double stepSize;
    virtual void operator()(std::shared_ptr<Particle> p, CallBack& move_callback) = 0;
    virtual bool accept(double dE) = 0;
};

class Translate : public Move{
    public:

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        std::vector< std::shared_ptr<Particle> > particles = {p};
        std::cout << "Before:" << std::endl;
        std::cout << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;
        p->translate();
        std::cout << "After:" << std::endl;
        std::cout << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;
        move_callback(particles);
    }

    bool accept(double dE){
        bool ret = false;
        (exp(-dE) > Random::get_random()) ? ret = false : ret = false;
        return ret;
    }
};

class Rotate : public Move{
    public:

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        std::vector< std::shared_ptr<Particle> > particles = {p};
        p->rotate();
        move_callback(particles);
    }

    bool accept(double dE){
        return exp(-dE);
    }
};
