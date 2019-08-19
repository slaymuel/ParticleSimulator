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
        p->translate();
        move_callback(particles);
    }

    bool accept(double dE){
        return exp(-dE);
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
