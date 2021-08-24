#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"
#include <algorithm>
#include <unordered_map>

using CallBack = std::function<void(std::vector< unsigned int >)>;



class Move{

    protected:
    int accepted, rejected, attempted;
    double stepSize;
    std::string id;
    State* s;
    CallBack move_callback;
    
    public:
    double weight;

    Move(double step, double w, State* state, CallBack move_callback) : stepSize(step), s(state), move_callback(move_callback), weight(w){
        this->accepted = 0;
        this->rejected = 0;
        this->attempted = 0;
    }

    virtual void operator()() = 0;
    virtual bool accept(double dE) = 0;
    virtual std::string dump() = 0;
};




class Translate : public Move{
    public:

    Translate(double step, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback){
        this->id = "Trans";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        std::shared_ptr<Particle> p = this->s->particles.random();
        //std::shared_ptr<Particle> p = std::static_pointer_cast<Particle>(argument);
        std::vector< unsigned int > particles = {p->index};
        //printf("Translating particle %lu\n", p->index);
        //std::cout << p->pos << std::endl;
        p->translate(this->stepSize);
        //printf("after move\n");
        //std::cout << p->pos << std::endl;
        //PBC


        this->move_callback(particles);
        this->attempted++;
    }

    bool accept(double dE){
        bool ret = false;
        //printf("dE trans %lf\n", dE);
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

    Rotate(double step, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback){
        this->id = "Rot";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        std::shared_ptr<Particle> p = this->s->particles.random();
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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};



class Swap : public Move{

    public:
    Swap(double w, State* s, CallBack move_callback) : Move(0.0, w, s, move_callback){
        this->id = "Swap";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }
    void operator()(){
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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};





class SingleSwap : public Move{

    public:
    SingleSwap(double w, State* s, CallBack move_callback) : Move(0.0, w, s, move_callback){
        //printf("\t\tSwap Move\n");
        this->id = "SingleSwap";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }
    
    void operator()(){

            
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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};





template <bool ADD>
class GrandCanonicalSingle : public Move{

    protected:
    double q;
    double d = 0.0;
    double nVolume;
    double pVolume;
    double cp;
    int pAtt = 0, nAtt = 0, pAcc = 0, nAcc = 0;
    int* att;
    int* acc;

    public:

    GrandCanonicalSingle(double chemPot, double donnan, double w, State* s, CallBack move_callback) : Move(0.0, w, s, move_callback),
                  d(donnan), cp(chemPot){
        if(ADD){
            this->id = "GCAddSingle";
        }
        else{
            this->id = "GCRemSingle";
        }
        printf("\t%s\n", this->id.c_str());
        constants::cp = chemPot;
        this->pVolume = this->s->geo->_d[0] * this->s->geo->_d[1] * (this->s->geo->_d[2] - 2.0 * this->s->particles.pModel.rf);
        this->nVolume = this->s->geo->_d[0] * this->s->geo->_d[1] * (this->s->geo->_d[2] - 2.0 * this->s->particles.nModel.rf);
        printf("\tCation accessible volume: %.3lf, Anion accessible volume: %.3lf\n", this->pVolume, this->nVolume);
        printf("\tChemical potential: %.3lf, Bias potential: %.3lf\n", this->cp, this->d);
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){
        //UNUSED(p);
        
        std::vector< unsigned int > particles;
        if(ADD){
            
            auto [ind, qt] = s->particles.add_random(s->geo->_dh);
            this->q = qt;
            particles.push_back(ind);
            this->move_callback(particles);
        }

        else{
            auto [ind, qt] = s->particles.remove_random();
            this->q = qt;
            if(ind != -1){
                particles.push_back(ind);
            }
            this->move_callback(particles);        
        }
        this->attempted++;
    }

    bool accept(double dE){
        double prob;

        // ADD
        if(ADD){
            //printf("dE add %lf\n", dE);
            //Cation
            if(this->q > 0.0){
                prob = this->pVolume / s->particles.cTot * std::exp(this->cp - this->d * this->q - dE); //N + 1 since s->particles.cTot is the new N + 1 state
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
            //printf("dE remove %lf\n", dE);
            //Cation
            if(this->q > 0.0){
                prob = (s->particles.cTot + 1) / this->pVolume * std::exp(this->d * this->q - this->cp - dE); //N since s->particles.cTot is the new N - 1 state
                this->acc = &this->pAcc;
                this->pAtt++;
            }
            //Anion
            else{
                prob = (s->particles.aTot + 1) / this->nVolume * std::exp(this->d * this->q - this->cp - dE);
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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << " +: " << (double) this->pAcc / this->pAtt * 100.0 << "%, " << this->pAtt << " (" << this->pAcc <<") ";
        ss << this->id << " -: " << (double) this->nAcc / this->nAtt * 100.0 << "%, " << this->nAtt << " (" << this->nAcc <<") ";

        return ss.str();
    }
};

















template <bool ADD>
class GrandCanonical : public Move{

    protected:
    unsigned int valency;
    double volume;
    double cp;
    double pVolume;
    double nVolume;

    public:

    GrandCanonical(double chemPot, double w, State* s, CallBack move_callback) : Move(0.0, w, s, move_callback), cp(chemPot){
        if(ADD){
            this->id = "GCAdd";
        }
        else{
            this->id = "GCRem";
        }
        
        constants::cp = chemPot;
        this->pVolume = this->s->geo->_d[0] * this->s->geo->_d[1] * (this->s->geo->_d[2] - 2.0 * this->s->particles.pModel.rf);
        this->nVolume = this->s->geo->_d[0] * this->s->geo->_d[1] * (this->s->geo->_d[2] - 2.0 * this->s->particles.nModel.rf);
        this->valency = s->particles.pModel.q;

        printf("\t%s\n", this->id.c_str());
        printf("\tCation accessible volume: %.3lf, Anion accessible volume: %.3lf\n", this->pVolume, this->nVolume);
        printf("\tChemical potential: %.3lf", this->cp);
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){
        //UNUSED(p);
        std::vector< unsigned int > particles;
        if(ADD){
            auto [ind, qt] = s->particles.add_random(s->geo->_dh, 1.0);
            particles.push_back(ind);
            for(int i = 0; i < this->valency; i++){
                auto [ind2, qt2] = s->particles.add_random(s->geo->_dh, -1.0);
                particles.push_back(ind2);
            }
        }
        else{
            std::vector< unsigned int > pt;
            /*int di;
            auto [ind, qt] = s->particles.remove_random(1.0);
            particles.push_back(ind);
            printf("removing: %i\n", ind);
            pt.push_back(ind);
            for(int i = 0; i < this->valency; i++){
                di = 0;
                auto [ind2, qt2] = s->particles.remove_random(-1.0);
                for(int j = 0; j < particles.size(); j++){
                    if(pt[j] <= ind2) di += 1;
                    if(pt[j] == ind2) di += 1;
                }
                printf("removing: %i\n", ind2);
                pt.push_back(ind2);
                particles.push_back(ind2 + di);
            }*/
            int rand = 0;
            if(s->particles.cTot > 0){
                do{
                    rand = Random::get_random(s->particles.tot);
                } while(s->particles[rand]->q != s->particles.pModel.q);
                particles.push_back(rand);
            }


            /*for(int i = 0; i < this->valency; i++){
                if(s->particles.aTot > this->valency){
                    do{
                        rand = Random::get_random(s->particles.tot);
                    } while(s->particles[rand]->q != s->particles.nModel.q && std::none_of(particles.begin(), particles.end(), [&](int val){return val==rand;}));
                    particles.push_back(rand);
                }
            }*/
            while(particles.size() < this->valency + 1){
                do{
                    rand = Random::get_random(s->particles.tot);
                } while(s->particles[rand]->q != s->particles.nModel.q);

                if(std::none_of(particles.begin(), particles.end(), [&](int val){return val==rand;})){
                    particles.push_back(rand);   
                }   
            }
            std::sort(particles.begin(), particles.end());
            std::reverse(particles.begin(), particles.end());
            for(auto p : particles){
                //printf("%i\n", p);
                s->particles.remove(p);
            }
            //printf("\n");
        }

        this->attempted++;
        this->move_callback(particles);  
    }

    bool accept(double dE){
        double prob;

        // ADD
        if(ADD){
            double pFac = this->pVolume / s->particles.cTot;
            double nFac = 1.0;
            for(int i = 0; i < this->valency; i++){
                nFac *= this->nVolume / (s->particles.aTot - i);
            }

            prob = pFac * nFac * std::exp((this->valency + 1) * this->cp - dE); 

        }
        // REMOVE
        else{
            double pFac = (s->particles.cTot + 1) / this->pVolume;
            double nFac = 1.0;
            for(int i = 0; i < this->valency; i++){
                nFac *= (s->particles.aTot + 1 - i) / this->nVolume ;
            }
            prob = pFac * nFac * std::exp(-(this->valency + 1.0) * this->cp - dE); 
        }

        if(prob >= Random::get_random()){
            this->accepted += 1;
            return true;
        }
        else{
            return false;
        }
    }

    std::string dump(){
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted/ this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted<<") ";

        return ss.str();
    }
};


















class VolumeMove: public Move{
    private:
    double _oldV;
    double pressure;
    double unit = 2.430527863808942e-10;

    public:
    VolumeMove(double step, double pressure, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback){
        this->id = "Vol";
        printf("\t%s\n", this->id.c_str());
        this->pressure = pressure * this->unit;
        printf("\tPressure: %lf\n", this->pressure);
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        _oldV = this->s->geo->volume;
        double lnV = std::log(this->s->geo->volume) + (Random::get_random() * 2.0 - 1.0) * this->stepSize;
        double V = std::exp(lnV);
        double L = std::cbrt(V);
        double RL = L / this->s->geo->_d[0];

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
        
        double prob = exp(-dE - this->pressure * (this->s->geo->volume - _oldV) +
                      (this->s->particles.tot + 1.0) * std::log(this->s->geo->volume / _oldV));

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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};




class ChargeTrans: public Move{

    public:

    ChargeTrans(double step, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback){
        this->id = "qTrans";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        int rand = 0;
        do{
            rand = Random::get_random(s->particles.tot);
        } while(s->particles[rand]->q < 0.0);
        
        std::vector< unsigned int > particles = {s->particles[rand]->index};
        //printf("Translating\n");
        this->s->particles[rand]->chargeTrans(this->stepSize);
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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};




class ChargeTranslate: public Move{

    public:

    ChargeTranslate(double step, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback){
        this->id = "qTranslate";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        int rand = 0;
        do{
            rand = Random::get_random(s->particles.tot);
        } while(s->particles[rand]->q < 0.0);
        
        std::vector< unsigned int > particles = {s->particles[rand]->index};
        this->s->particles[rand]->chargeTranslate(this->stepSize);

        this->s->particles[rand]->qDisp = this->s->geo->displacement(this->s->particles[rand]->pos, this->s->particles[rand]->com);
        this->s->particles[rand]->b = this->s->particles[rand]->qDisp.norm();

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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};



class ChargeTransRand: public Move{

    public:

    ChargeTransRand(double step, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback){
        this->id = "qTransRand";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        //printf("Move\n");
        std::vector< unsigned int > particles;
        int rand = 0;
        if(s->particles.cTot > 0){
            do{
                rand = Random::get_random(s->particles.tot);
            } while(s->particles[rand]->q < 0.0);

            this->s->particles[rand]->chargeTransRand();
            //this->s->particles[rand]->qDisp = this->s->geo->displacement(this->s->particles[rand]->pos, this->s->particles[rand]->com);
            //this->s->particles[rand]->b = this->s->particles[rand]->qDisp.norm();
            
            particles.push_back(s->particles[rand]->index);
        }
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
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};


class Cluster : public Move{
    private:
    double minDist;
    bool found = false;
    std::unordered_map<unsigned int, unsigned int> att;
    std::unordered_map<unsigned int, unsigned int> acc;
    unsigned int pNum = 0;
    std::shared_ptr<Particle> p;

    public:

    Cluster(double step, double md, double w, State* s, CallBack move_callback) : Move(step, w, s, move_callback), minDist(md){
        this->id = "Clus";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tMax distance in cluster: %lf\n", this->minDist);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        this->p = this->s->particles.random();
        std::vector<unsigned int> indices;
        Eigen::Vector3d disp;

        for(auto i : this->s->particles.particles){
            if(i->index == this->p->index) continue;

            if(this->s->geo->distance(i->pos, this->p->pos) <= this->minDist){
                indices.push_back(i->index);
            }
        }
        //Choose random particle
        //If particle is in cluster
        //Move cluster

        if(!indices.empty()){
            indices.push_back(this->p->index);
            this->pNum = indices.size();

            disp = Random::get_norm_vector();
            disp *= this->stepSize;
            this->s->particles.translate(indices, disp);
            this->move_callback(indices);
            this->attempted++;
            found = true;

            if(att.count(this->pNum) > 0){
                att[this->pNum]++;
            }
            else{
                att.insert(std::pair<int, int>(this->pNum, 1));
                acc.insert(std::pair<int, int>(this->pNum, 0));
            }
        }
        else found = false;
    }

    bool accept(double dE){
        int count = 1;

        if(found){

            for(auto i : this->s->particles.particles){
                if(i->index == this->p->index) continue;
                if(this->s->geo->distance(i->pos, this->p->pos) <= this->minDist){
                    count++;
                }
            }
            if(count != this->pNum) return false;
            
            if(exp(-dE) >= Random::get_random() || dE < 0.0){
                acc[this->pNum]++;
                this->accepted++;
                return true;
            } 
            else{
                this->rejected++;
                return false;
            }
        }
        return false;
    }

    std::string dump(){
        std::ostringstream ss;
        ss.precision(1);
        ss << std::fixed;
        ss << "\t" << this->id << " ";
        for(auto const& [key, val] : this->att){
            ss << key << ": " << (double)this->acc[key] / val * 100.0 << "%, " << val << "   ";
        } 
        ss <<" Total: " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
        return ss.str();
    }
};


class WidomInsertion : public Move{

    protected:
    double q;
    double cp = 0.0;
    int samples = 0;
    public:

    WidomInsertion(double w, State* s, CallBack move_callback) : Move(0.0, w, s, move_callback) {

        this->id = "WIns";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){      
        auto [ind, qt] = s->particles.add_random(s->geo->_dh);
        this->q = qt;
        std::vector< unsigned int > particles{ind};
        this->move_callback(particles);

        this->attempted++;
    }

    bool accept(double dE){
        if(this->samples > 10000){
            this->samples = 0;
            this->cp = 0.0;
        }
        
        this->cp += std::exp(-dE);
        this->samples++;

        return false;
    }

    std::string dump(){
        std::ostringstream ss;
        ss.precision(5);
        ss << std::fixed;
        ss << "\t" << this->id << " cp: " << -std::log(this->cp / this->samples) <<" " << this->cp / this->samples << " " << this->cp << " samples: " << this->samples << " attempted: " << this->attempted;
        return ss.str();
    }
};

class WidomDeletion : public Move{

    protected:
    double q;
    double cp = 0.0;
    int samples = 0;
    public:

    WidomDeletion(double w, State* s, CallBack move_callback) : Move(0.0, w, s, move_callback) {

        this->id = "WDel";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){      
        auto [ind, qt] = s->particles.remove_random();
        this->q = qt;
        std::vector< unsigned int > particles{ind};
        this->move_callback(particles);

        this->attempted++;
    }

    bool accept(double dE){
        if(this->samples > 10000){
            this->samples = 0;
            this->cp = 0.0;
        }
        //printf("dE %lf\n", dE);
        this->cp += std::exp(dE);
        this->samples++;

        return false;
    }

    std::string dump(){
        std::ostringstream ss;
        ss.precision(5);
        ss << std::fixed;
        ss << "\t" << this->id << " cp: " << -std::log(this->cp / this->samples) <<" " << this->cp / this->samples << " " << this->cp << " samples: " << this->samples << " attempted: " << this->attempted;
        return ss.str();
    }
};