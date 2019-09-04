#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"

using CallBack = std::function<void(std::vector< int >)>;

class Move{
    public:
    State* s;
    int accepted, rejected, attempted;
    double stepSize;

    static int totalMoves;

    Move(double step) : stepSize(step){
        accepted = 0;
        rejected = 0;
        attempted = 0;
    }
    virtual void operator()(std::shared_ptr<void> argument, CallBack& move_callback) = 0;
    virtual bool accept(double dE) = 0;
};



template <bool HW>
class Translate : public Move{
    public:

    Translate(double step) : Move(step){}


    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< int > particles = {p->index};
        //printf("Translating\n");
        p->translate(this->stepSize);
        //PBC

        if(HW){
            //Move image particle
        }

        move_callback(particles);
        totalMoves++;
        attempted++;
    }

    bool accept(double dE){
        bool ret = false;
        if(exp(-dE) > Random::get_random()){
            ret = true;
            accepted++;
            //printf("Accepted\n");
         } 
         else{
            ret = false;
            rejected++;
            //printf("Rejected\n");
         }
        return ret;
    }
};



class Rotate : public Move{
    public:

    Rotate(double step) : Move(step){}


    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< int > particles = {p->index};

        p->rotate(this->stepSize);
        //PBC
        
        move_callback(particles);
        totalMoves++;
        attempted++;
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
class GrandCanonicalAdd : public Move{
    private:
    double cp, d;

    public:
    GrandCanonicalAdd(double chemPot, double donnan = 0) :  Move(0.0), cp(chemPot), d(donnan){}

    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        //std::shared_ptr<State> s = std::static_pointer_cast<State>(argument);
        //printf("Adding\n");
        s->particles.add(s->geo->random_pos(), 2.5, 1, 0, "add"); // add particle
        std::vector< int > particles{s->particles.tot - 1};
        //printf("Move: Added particle %i\n", s->particles.tot - 1);
        move_callback(particles);

        if(HW){
            //p->add() // add image
        }

        totalMoves++;
        attempted++;
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
            //printf("Accepted\n");
         } 
         else{
             ret = false;
             rejected++;
             //printf("Rejected\n");
         }
        return ret;
    }
};

template <bool HW>
class GrandCanonicalRemove : public Move{
    private:
    double cp, d;

    public:
    GrandCanonicalRemove(double chemPot, double donnan = 0) :  Move(0.0), cp(chemPot), d(donnan){}

    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< int > particles = {p->index};
        //printf("\nremoving particle %i from current\n\n", p->index);
        //printf("Removing\n");
        s->particles.remove(p->index); // add particle

        if(HW){
            // remove image
        }

        move_callback(particles);

        totalMoves++;
        attempted++;
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
            //printf("Accepted\n");
         } 
         else{
             ret = false;
             rejected++;
             //printf("Rejected\n");
         }
        return ret;
    }
};


int Move::totalMoves = 0;
