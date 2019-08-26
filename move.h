#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>

using CallBack = std::function<void(std::vector< std::shared_ptr<Particle> >)>;

class Move{
    public:

    int accepted, rejected;
    double stepSize;
    static int totalMoves;
    ;

    Move(double step) : stepSize(step){
        accepted = 0;
        rejected = 0;
    }
    virtual void operator()(std::shared_ptr<void> argument, CallBack& move_callback) = 0;
    virtual bool accept(double dE) = 0;
};


class Translate : public Move{
    public:

    Translate(double step) : Move(step){}


    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< std::shared_ptr<Particle> > particles = {p};

        p->translate(this->stepSize);

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


    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< std::shared_ptr<Particle> > particles = {p};

        p->rotate(this->stepSize);

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


template <bool HW>
class GrandCanonical : public Move{
    private:
    double cp, d;

    public:
    GrandCanonical(double chemPot, double donnan = 0) : cp(chemPot), d(donnan){}

    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        std::shared_ptr<Particles> p = std::static_pointer_cast<Particles>(argument);

        //p->add(); // add particle
        if(HW){
            //p->add() // add image
        }
    }


    bool accept(double dE){
        //add anion
        //prob =  volumeN / particles.numOfAnions * std::exp(cp - d * q - dE);
        //Add cation
        //prob = volumeP / particles.numOfCations * std::exp(cp - d * q - dE);
        //delete anion
        //prob =  particles.numOfAnions / volumeN * std::exp(d * particles[r].q - cp + dE);
        //delete cation
        //prob = particles.numOfCations / volumeP * std::exp(d * particles[r].q - cp + dE);

    }
};



int Move::totalMoves = 0;
