#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"

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
    State* s;
    GrandCanonical(double chemPot, double donnan = 0) :  Move(0.0), cp(chemPot), d(donnan){}

    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        //std::shared_ptr<State> s = std::static_pointer_cast<State>(argument);
        
        s->particles.add(s->geo->random_pos(), 2.5, 1, 0, "add"); // add particle
        std::vector< std::shared_ptr<Particle> > particles = {s->particles.particles[s->particles.tot - 1]};
        move_callback(particles);

        if(HW){
            //p->add() // add image
        }

        totalMoves++;
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



int Move::totalMoves = 0;
