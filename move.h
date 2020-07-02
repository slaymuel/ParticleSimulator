#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"
#include <algorithm>

using CallBack = std::function<void(std::vector< unsigned int >)>;



class Move{

    protected:
    int accepted, rejected, attempted;
    double stepSize;
    std::string id;
    CallBack move_callback;

    public:
    double weight;

    Move(double step, double w, CallBack move_callback) : stepSize(step), move_callback(move_callback), weight(w){
        this->accepted = 0;
        this->rejected = 0;
        this->attempted = 0;
    }

    virtual void operator()(std::shared_ptr<Particle> p) = 0;
    virtual bool accept(double dE) = 0;
    virtual std::string dump() = 0;
};




class Translate : public Move{
    public:

    Translate(double step, double w, CallBack move_callback) : Move(step, w, move_callback){
        this->id = "Trans";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
    }


    void operator()(std::shared_ptr<Particle> p){
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< unsigned int > particles = {p->index};
        //printf("Translating\n");
        p->translate(this->stepSize);
        //PBC


        this->move_callback(particles);
        this->attempted++;
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }
    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};



class Rotate : public Move{
    public:

    Rotate(double step, double w, CallBack move_callback) : Move(step, w, move_callback){
        this->id = "Rot";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
    }


    void operator()(std::shared_ptr<Particle> p){
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< unsigned int > particles = {p->index};

        p->rotate(this->stepSize);
        this->move_callback(particles);
        attempted++;
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }
    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};



class Swap : public Move{
    private:
    State* s;

    public:
    Swap(State* state, double w, CallBack move_callback) : Move(0.0, w, move_callback), s(state){
        this->id = "Swap";
        printf("\t%s\n", this->id.c_str());
    }
    void operator()(std::shared_ptr<Particle> p){
        int rand = Random::get_random(s->particles.tot), rand2;

        do{
            rand2 = Random::get_random(s->particles.tot);
        } while(this->s->particles[rand]->q == this->s->particles[rand2]->q);

        /*std::swap(this->s->particles.particles[rand]->q, this->s->particles.particles[rand2]->q);
        std::swap(this->s->particles.particles[rand]->name, this->s->particles.particles[rand2]->name);
        std::swap(this->s->particles.particles[rand]->r, this->s->particles.particles[rand2]->r);
        std::swap(this->s->particles.particles[rand]->rf, this->s->particles.particles[rand2]->rf);
        std::swap(this->s->particles.particles[rand]->b, this->s->particles.particles[rand2]->b);
        std::swap(this->s->particles.particles[rand]->qDisp, this->s->particles.particles[rand2]->qDisp);*/

        std::swap(this->s->particles.particles[rand]->pos, this->s->particles.particles[rand2]->pos);
        std::swap(this->s->particles.particles[rand]->com, this->s->particles.particles[rand2]->com);
        std::vector< unsigned int > particles = {static_cast<unsigned int>(rand), static_cast<unsigned int>(rand2)};


        this->move_callback(particles);
        attempted++;

        //printf("Remove: pAtt: %i, nAtt: %i\n", this->pAtt, this->nAtt);
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }

    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};





class SingleSwap : public Move{
    private:
    State* s;

    public:
    SingleSwap(State* state, double w, CallBack move_callback) : Move(0.0, w, move_callback), s(state){
        //printf("\t\tSwap Move\n");
        this->id = "SingleSwap";
        printf("\t%s\n", this->id.c_str());
    }
    
    void operator()(std::shared_ptr<Particle> p){

            
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        //std::vector< unsigned int > particles = {p->index};
        int rand = Random::get_random(s->particles.tot);
        //int rand2;

        //If cation
        if(this->s->particles[rand]->q > 0){
            this->s->particles[rand]->q = this->s->particles.nModel.q;
            this->s->particles[rand]->b = this->s->particles.nModel.b;
            this->s->particles[rand]->r = this->s->particles.nModel.r;
            this->s->particles[rand]->rf = this->s->particles.nModel.rf;

            //Set qDisp
            Eigen::Vector3d v = Random::get_vector();
            this->s->particles[rand]->qDisp = v;
            this->s->particles[rand]->qDisp = this->s->particles[rand]->qDisp.normalized() * this->s->particles[rand]->b;
            this->s->particles[rand]->pos = this->s->particles[rand]->com + this->s->particles[rand]->qDisp;

            this->s->particles[rand]->name = "Cl";
            this->s->particles.aTot++;
            this->s->particles.cTot--;
        }

        //anion
        else{
            //flip charge and change name
            this->s->particles[rand]->q = this->s->particles.pModel.q;
            this->s->particles[rand]->b = this->s->particles.pModel.b;
            this->s->particles[rand]->r = this->s->particles.pModel.r;
            this->s->particles[rand]->rf = this->s->particles.pModel.rf;

            //Set qDisp
            Eigen::Vector3d v = Random::get_vector();
            this->s->particles[rand]->qDisp = v;
            this->s->particles[rand]->qDisp = this->s->particles[rand]->qDisp.normalized() * this->s->particles[rand]->b;
            this->s->particles[rand]->pos = this->s->particles[rand]->com + this->s->particles[rand]->qDisp;

            this->s->particles[rand]->name = "Na";
            this->s->particles.cTot++;
            this->s->particles.aTot--;      
        }


        std::vector< unsigned int > particles = {static_cast<unsigned int>(rand)};


        this->move_callback(particles);
        attempted++;
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }

    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};





template <bool ADD>
class GrandCanonical : public Move{

    protected:
    State* s;
    double q;
    double d = 0.0;
    double nVolume;
    double pVolume;
    double cp;
    int pAtt = 0, nAtt = 0, pAcc = 0, nAcc = 0;
    int* att;
    int* acc;

    public:

    GrandCanonical(double chemPot, double donnan, State* state, double w, CallBack move_callback) : Move(0.0, w, move_callback), s(state), d(donnan){
        if(ADD){
            this->id = "GCAdd";
        }
        else{
            this->id = "GCRem";
        }
        printf("\t%s\n", this->id.c_str());
        constants::cp = chemPot;
        this->cp = chemPot;
        this->pVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.pModel.rf);
        this->nVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.nModel.rf);
        printf("\t\tCation accessible volume: %.3lf, Anion accessible volume: %.3lf\n", this->pVolume, this->nVolume);
        printf("\t\tChemical potential: %.3lf, Bias potential: %.3lf\n", constants::cp, this->d);
    }

    void operator()(std::shared_ptr<Particle> p){
        UNUSED(p);
        if(ADD){
            auto [ind, qt] = s->particles.add_random(s->geo->_dh);
            this->q = qt;
            std::vector< unsigned int > particles{ind};
            this->move_callback(particles);
        }
        else{
            auto [ind, qt] = s->particles.remove_random();
            this->q = qt;
            std::vector< unsigned int > particles{ind};
            this->move_callback(particles);        
        }
        this->attempted++;
    }

    bool accept(double dE){
        double prob;

        // ADD
        if(ADD){
            //Cation
            if(this->q > 0.0){
                prob = this->pVolume / s->particles.cTot * std::exp(this->cp - this->d * this->q - dE);
                this->acc = &this->pAcc;
                this->pAtt++;
            }
            //Anion
            else{
                prob = this->nVolume / s->particles.aTot * std::exp(this->cp - this->d * this->q - dE);
                this->acc = &this->nAcc;
                this->nAtt++;
            } 
        }
        // REMOVE
        else{
            //Cation
            if(this->q > 0.0){
                prob = (s->particles.cTot + 1.0) / this->pVolume * std::exp(this->d * this->q - this->cp - dE);
                this->acc = &this->pAcc;
                this->pAtt++;
            }
            //Anion
            else{
                prob = (s->particles.aTot + 1.0) / this->nVolume * std::exp(this->d * this->q - this->cp - dE);
                this->acc = &this->nAcc;
                this->nAtt++;
            }
        }

        if(prob >= Random::get_random()){
            *(this->acc) += 1;
            return true;
        }
        else{
            return false;
        }
    }

    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << " +: " << (double) this->pAcc / this->pAtt * 100.0 << "%, " << this->pAtt << " (" << this->pAcc <<") ";
        s << this->id << " -: " << (double) this->nAcc / this->nAtt * 100.0 << "%, " << this->nAtt << " (" << this->nAcc <<") ";

        return s.str();
    }
};





class VolumeMove: public Move{
    private:
    State* s;
    double _oldV;
    double pressure;
    double unit = 2.430527863808942e-10;

    public:
    VolumeMove(State* s, double step, double pressure, double w, CallBack move_callback) : s(s), Move(step, w, move_callback){
        this->id = "Vol";
        printf("\t%s\n", this->id.c_str());
        this->pressure = pressure * this->unit;
        printf("\tPressure: %lf\n", this->pressure);
        printf("\tStepsize: %lf\n", step);
    }


    void operator()(std::shared_ptr<Particle> p){
        _oldV = this->s->geo->volume;
        //double vMax = 0.00025;
        double lnV = std::log(this->s->geo->volume) + (Random::get_random() * 2.0 - 1.0) * this->stepSize;
        double V = std::exp(lnV);
        //V = this->s->geo->volume + (Random::get_random() * 2.0 - 1.0) * this->stepSize;
        //printf("Changing volume by: %lf\n", V - this->s->geo->volume);
        double L = std::cbrt(V);
        double RL = L / this->s->geo->_d[0];
        //double oldL = this->s->geo->_d[0];

        std::vector<double> LV = {L, L, L};
        std::vector<double> LVh = {L / 2.0, L / 2.0, L / 2.0};

        this->s->geo->_d = LV;
        this->s->geo->_dh = LVh;
        this->s->geo->d = LV;
        this->s->geo->dh = LVh;
        this->s->geo->volume = V;

        std::vector< unsigned int > particles;

        for(unsigned int i = 0; i < this->s->particles.tot; i++){
            this->s->particles[i]->com *= RL;
            this->s->particles[i]->pos = this->s->particles[i]->com + this->s->particles[i]->qDisp;
            particles.push_back(s->particles[i]->index);
        }

        this->move_callback(particles);
        this->attempted++;
    }

    bool accept(double dE){
        bool ret = false;
        
        double prob = exp(-dE - this->pressure * (this->s->geo->volume - _oldV) + //  0.0000243     0.005      0.00383374 0.000024305278638
                      (this->s->particles.tot + 1) * std::log(this->s->geo->volume / _oldV));

        if(prob >= Random::get_random()){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }

    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};




class ChargeTrans: public Move{
    private:
    State* s;
    public:

    ChargeTrans(State* s, double step, double w, CallBack move_callback) : s(s), Move(step, w, move_callback){
        this->id = "qTrans";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
    }


    void operator()(std::shared_ptr<Particle> p){
        int rand = 0;
        do{
            rand = Random::get_random(s->particles.tot);
        } while(s->particles[rand]->q < 0.0);
        
        std::vector< unsigned int > particles = {s->particles[rand]->index};
        //printf("Translating\n");
        s->particles[rand]->chargeTrans(this->stepSize);
        this->move_callback(particles);
        this->attempted++;
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }

    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};




class ChargeTransRand: public Move{
    private:
    State* s;
    public:

    ChargeTransRand(State* s, double step, double w, CallBack move_callback) : s(s), Move(step, w, move_callback){
        this->id = "qTransRand";
        printf("\t%s\n", this->id.c_str());
    }


    void operator()(std::shared_ptr<Particle> p){
        //printf("Move\n");
        int rand = 0;
        do{
            rand = Random::get_random(s->particles.tot);
        } while(s->particles[rand]->q < 0.0);
        
        std::vector< unsigned int > particles = {s->particles[rand]->index};
        s->particles[rand]->chargeTransRand();
        this->move_callback(particles);
        this->attempted++;
    }

    bool accept(double dE){
        bool ret = false;

        if(exp(-dE) >= Random::get_random() || dE < 0.0){
            ret = true;
            this->accepted++;
         } 
         else{
            ret = false;
            this->rejected++;
         }

        return ret;
    }

    std::string dump(){
        std::ostringstream s;
        s.precision(1);
        s << std::fixed;
        s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return s.str();
    }
};