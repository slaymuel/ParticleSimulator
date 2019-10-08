#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"

using CallBack = std::function<void(std::vector< unsigned int >)>;

class Move{
    public:
    int accepted, rejected, attempted;
    double stepSize;
    double weight;
    std::string id;

    static int totalMoves;

    Move(double step, double w) : stepSize(step), weight(w){
        accepted = 0;
        rejected = 0;
        attempted = 0;
    }
    virtual void operator()(std::shared_ptr<Particle> p, CallBack& move_callback) = 0;
    virtual bool accept(double dE) = 0;
};




class Translate : public Move{
    public:

    Translate(double step, double w) : Move(step, w){
        this->id = "Trans";
    }


    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< unsigned int > particles = {p->index};
        //printf("Translating\n");
        p->translate(this->stepSize);
        //PBC


        move_callback(particles);
        totalMoves++;
        this->attempted++;
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
            //printf("Accepted\n");
         } 
         else{
            ret = false;
            this->rejected++;
            //printf("Rejected\n");
         }

        return ret;
    }
};



class Rotate : public Move{
    public:

    Rotate(double step, double w) : Move(step, w){
        this->id = "Rot";
    }


    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< unsigned int > particles = {p->index};

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
    double cp = 0.0, d = 0.0;
    double nVolume;
    double pVolume;
    int pAtt = 0, nAtt = 0, pAcc = 0, nAcc = 0;

    public:

    GrandCanonical(double chemPot, double donnan, State* state, double w) : Move(0.0, w), s(state) ,cp(chemPot), d(donnan){
        this->pVolume = state->geo->d[0] * state->geo->d[1] * (state->geo->dh[2] - 2.0 * state->particles.pModel.rf);
        this->nVolume = state->geo->d[0] * state->geo->d[1] * (state->geo->dh[2] - 2.0 * state->particles.nModel.rf);
        printf("\t\tCation accessible volume: %.3lf, Anion accessible volume: %.3lf\n", this->pVolume, this->nVolume);
        printf("\t\tChemical potential: %.3lf, Bias potential: %.3lf\n", this->cp, this->d);
    }

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback) = 0;
    bool accept(double dE) = 0;

    template<bool ADD>
    bool accept_imp(double dE){

        double prob;

        if(ADD){
            if(this->q > 0){
                prob = this->pVolume / s->particles.cTot * std::exp(this->cp - this->d * this->q - dE);;
                //prob = std::exp(-dE);
                //prob = volumeP / particles.numOfCations * std::exp(cp - d * q - dE);
            }
            else{
                prob = this->nVolume / s->particles.aTot * std::exp(this->cp - this->d * this->q - dE);
                //prob = std::exp(-dE);
                //prob =  volumeN / particles.numOfAnions * std::exp(cp - d * q - dE);
            }            
        }
        else{
            if(this->q > 0){
                prob = s->particles.cTot / this->pVolume * std::exp(this->d * this->q - this->cp - dE);
                //particles.numOfCations / volumeP * std::exp(d * particles[r].q - cp + dE);
            }
            else{
                prob = s->particles.aTot / this->nVolume * std::exp(this->d * this->q - this->cp - dE);
                //particles.numOfAnions / volumeN * std::exp(d * particles[r].q - cp + dE);
            }
        }
        //printf("%lf\n", dE);
        if(prob >= Random::get_random()){
            return true;
        }
        else{
            return false;
        }
    }
};


class GrandCanonicalAdd : public GrandCanonical{

    public:
    GrandCanonicalAdd(double chemPot, double donnan, State* state, double w) : GrandCanonical(chemPot, donnan, state, w){
        this->id = "GCAdd";
    }

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){
        UNUSED(p);
        //std::shared_ptr<State> s = std::static_pointer_cast<State>(argument);
        double rand = Random::get_random();
        if(rand < 0.5){//if(rand < s->particles.cTot / s->particles.tot){
            s->particles.add(s->geo->random_pos(), s->particles.pModel.r, s->particles.pModel.rf, s->particles.pModel.q, 0, "Na");
            this->q = s->particles.pModel.q;
            this->pAtt++;
        }
        else{
            s->particles.add(s->geo->random_pos(), s->particles.nModel.r, s->particles.nModel.rf, s->particles.nModel.q, 0, "Cl");
            this->q = s->particles.nModel.q;
            this->nAtt++;
        }

        std::vector< unsigned int > particles{s->particles.tot - 1};

        move_callback(particles);
        totalMoves++;
        attempted++;
        //printf("Add: pAtt: %i, nAtt: %i\n", this->pAtt, this->nAtt);
    }


    bool accept(double dE){
        if(accept_imp<true>(dE)){
            this->accepted++;
            return true;
        }
        else{
            this->rejected++;
            return false;
        }
    }
};


class GrandCanonicalRemove : public GrandCanonical{
    public:
    GrandCanonicalRemove(double chemPot, double donnan, State* state, double w) : GrandCanonical(chemPot, donnan, state, w){
        this->id = "GCRem";
    }

    void operator()(std::shared_ptr<Particle> p, CallBack& move_callback){

            
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        //std::vector< unsigned int > particles = {p->index};
        double rand = Random::get_random();
        int rand2 = 0;
        if(rand < 0.5){
            if(s->particles.cTot > 0){
                do{
                    rand2 = Random::get_random(s->particles.tot);
                } while(s->particles[rand2]->q != s->particles.pModel.q);
                pAtt++;
            }
        }
        else{
            if(s->particles.aTot > 0){
                do{
                    rand2 = Random::get_random(s->particles.tot);
                } while(s->particles[rand2]->q != s->particles.nModel.q);
                nAtt++;
            }
        }

        p = s->particles[rand2];
        std::vector< unsigned int > particles = {p->index};
        this->q = p->q;

        s->particles.remove(p->index); // remove particle

        move_callback(particles);
        totalMoves++;
        attempted++;

        //printf("Remove: pAtt: %i, nAtt: %i\n", this->pAtt, this->nAtt);
    }

    bool accept(double dE){
        if(accept_imp<false>(dE)){
            this->accepted++;
            return true;
        }
        else{
            this->rejected++;
            return false;
        }
    }
};


int Move::totalMoves = 0;
