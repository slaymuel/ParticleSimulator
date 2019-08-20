#pragma once
#include "particle.h"
#include <math.h>

using CallBack = std::function<void(std::vector< std::shared_ptr<Particle> >)>;

class Move{
    public:

    int accepted, rejected;
    double stepSize;
    static int totalMoves;

    Move(double step) : stepSize(step){
        accepted = 0;
        rejected = 0;
    }
    virtual void operator()(std::shared_ptr<Particle> p, CallBack& move_callback) = 0;
    virtual bool accept(double dE) = 0;
};

class Translate : public Move{
    public:

    Translate(double step) : Move(step){}

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        std::vector< std::shared_ptr<Particle> > particles = {p};
        std::cout << "Before:" << std::endl;
        std::cout << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;

        p->translate(this->stepSize);

        std::cout << "After:" << std::endl;
        std::cout << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;
        move_callback(particles);
        totalMoves++;
    }

    bool accept(double dE){
        bool ret = false;
        if(exp(-dE) > Random::get_random()){
            ret = true;
            accepted++;
         } 
         else{
             ret = false;
             rejected++;
         }
        return ret;
    }
};

class Rotate : public Move{
    public:
    
    Rotate(double step) : Move(step){}

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        std::vector< std::shared_ptr<Particle> > particles = {p};

        p->rotate(this->stepSize);

        move_callback(particles);
    }

    bool accept(double dE){
        return exp(-dE);
    }
};

int Move::totalMoves = 0;
