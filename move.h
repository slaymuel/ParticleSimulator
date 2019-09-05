#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"

using CallBack = std::function<void(std::vector< int >)>;

class Move{
    public:
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

class GrandCanonical : public Move{

    protected:
    State* s;
    double q;
    double cp = 10.0, d = 0.0;
    
    public:

    GrandCanonical(double chemPot, double donnan, State* state) : Move(0.0), s(state) ,cp(chemPot), d(donnan){}
    virtual void operator()(std::shared_ptr<void> argument, CallBack& move_callback){};
    virtual bool accept(double dE){return false;};

    template<bool ADD>
    bool accept_imp(double dE){
        //delete anion

        double prob;
        
        if(ADD){
            if(this->q > 0){
                prob = s->geo->volume / s->particles.cTot * std::exp(cp - d * q - dE);;
                //prob = std::exp(-dE);
                //prob = volumeP / particles.numOfCations * std::exp(cp - d * q - dE);
            }
            else{
                prob = s->geo->volume / s->particles.aTot * std::exp(cp - d * q - dE);
                //prob = std::exp(-dE);
                //prob =  volumeN / particles.numOfAnions * std::exp(cp - d * q - dE);
            }            
        }
        else{
            if(this->q > 0){
                prob = s->particles.cTot / s->geo->volume * std::exp(d * this->q - cp - dE);
                //particles.numOfCations / volumeP * std::exp(d * particles[r].q - cp + dE);
            }
            else{
                prob =  s->particles.aTot / s->geo->volume * std::exp(d * this->q - cp - dE);
                //particles.numOfAnions / volumeN * std::exp(d * particles[r].q - cp + dE);
            }
        }
        //printf("%lf\n", dE);
        return prob > Random::get_random();
    }
};

template <bool HW>
class GrandCanonicalAdd : public GrandCanonical{

    public:
    GrandCanonicalAdd(double chemPot, double donnan, State* state) : GrandCanonical(chemPot, donnan, state){}

    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        //std::shared_ptr<State> s = std::static_pointer_cast<State>(argument);
        //printf("Adding\n");
        this->q = (Random::get_random() > 0.5) ? 1.0 : -1.0;
        s->particles.add(s->geo->random_pos(), 2.5, this->q, 0, (this->q > 0) ? "Na" : "Cl"); // add particle
        std::vector< int > particles{s->particles.tot - 1};

        //printf("Move: Added particle %i\n", s->particles.tot - 1);

        if(HW){
            //p->add() // add image
        }

        move_callback(particles);
        totalMoves++;
        attempted++;
    }


    bool accept(double dE){
        if(accept_imp<true>(dE)){
            accepted++;
            return true;
        }
        else{
            rejected++;
            return false;
        }
        
    }
};

template <bool HW>
class GrandCanonicalRemove : public GrandCanonical{
    public:
    GrandCanonicalRemove(double chemPot, double donnan, State* state) : GrandCanonical(chemPot, donnan, state){}

    void operator()(std::shared_ptr<void> argument, CallBack& move_callback){
        if(s->particles.tot > 0){
            std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
            std::vector< int > particles = {p->index};
            this->q = p->q;
            //printf("\nremoving particle %i from current\n\n", p->index);
            //printf("Removing\n");
            s->particles.remove(p->index); // remove particle

            if(HW){
                // remove image
            }

            move_callback(particles);
            totalMoves++;
            attempted++;
        }
        else{
            printf("All particles have been removed!\n");
            exit(1);
        }
    }

    bool accept(double dE){
        if(accept_imp<false>(dE)){
            accepted++;
            return true;
        }
        else{
            rejected++;
            return false;
        }
    }
};


int Move::totalMoves = 0;
