#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"
#include <algorithm>
#include <unordered_map>

namespace Simulator{

using CallBack = std::function<void(std::vector< unsigned int >)>;

enum class MoveTypes{
    Translate,
    GCSingleAdd,
    GCSingleRemove,
    Rotate,
    Swap,
    SingleSwap,
    VolumeMove,
    ChargeTrans,
    ChargeTransRand,
    Cluster,
    WidomInsertion,
    WidomDeletion,
    GCAdd,
    GCRemove,
    ChargeTranslate
};




class Move{

    protected:
    mutable int accepted, rejected, attempted;
    double stepSize;
    std::string id;

    public:
    double weight;
    static State* state;

    Move(double step, double w) : stepSize(step), weight(w){
        this->accepted = 0;
        this->rejected = 0;
        this->attempted = 0;
    }

    virtual void operator()() = 0;
    virtual bool accept(double dE) const = 0;
    virtual std::string dump() = 0;

    static Move* create_move(MoveTypes moveType, std::vector<double> args);
};

std::string Move::dump(){
    std::ostringstream s;
    s.precision(1);
    s << std::fixed;
    s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
    return s.str();
}

State* Move::state = nullptr;


class Translate : public Move{
    public:

    Translate(double step, double w) : Move(step, w){
        this->id = "Trans";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        auto p = state->particles.random();

        p->translate(stepSize);
        this->attempted++;

        state->move_callback({p->index});
    }

    bool accept(double dE) const{
        auto  ret = false;

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
        return Move::dump();
    }
};



class Rotate : public Move{
    public:

    Rotate(double step, double w) : Move(step, w){
        this->id = "Rot";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        auto p = state->particles.random();

        p->rotate(this->stepSize);
        attempted++;

        state->move_callback({p->index});
    }

    bool accept(double dE) const{
        auto ret = false;

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
        return Move::dump();
    }
};



class Swap : public Move{

    public:
    Swap(double w) : Move(0.0, w){
        this->id = "Swap";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }
    void operator()(){
        int rand2 = -1;
        auto rand = Random::get_random(state->particles.tot);

        do{
            rand2 = Random::get_random(state->particles.tot);
        } while(state->particles[rand]->q == state->particles[rand2]->q);

        std::swap(state->particles.particles[rand]->pos, state->particles.particles[rand2]->pos);
        std::swap(state->particles.particles[rand]->com, state->particles.particles[rand2]->com);

        attempted++;
        state->move_callback({static_cast<unsigned int>(rand), static_cast<unsigned int>(rand2)});
    }

    bool accept(double dE) const{
        auto ret = false;

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
        return Move::dump();
    }
};





class SingleSwap : public Move{

    public:
    SingleSwap(double w) : Move(0.0, w){
        //printf("\t\tSwap Move\n");
        this->id = "SingleSwap";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }
    
    void operator()(){
        auto rand = Random::get_random(state->particles.tot);

        //If cation
        if(state->particles[rand]->q > 0){
            state->particles[rand]->q = state->particles.nModel.q;
            state->particles[rand]->b = state->particles.nModel.b;
            state->particles[rand]->r = state->particles.nModel.r;
            state->particles[rand]->rf = state->particles.nModel.rf;

            //Set qDisp
            Eigen::Vector3d v = Random::get_vector();
            state->particles[rand]->qDisp = v;
            state->particles[rand]->qDisp = state->particles[rand]->qDisp.normalized() * state->particles[rand]->b;
            state->particles[rand]->pos = state->particles[rand]->com + state->particles[rand]->qDisp;

            state->particles[rand]->name = "Cl";
            state->particles.aTot++;
            state->particles.cTot--;
        }

        //anion
        else{
            //flip charge and change name
            state->particles[rand]->q = state->particles.pModel.q;
            state->particles[rand]->b = state->particles.pModel.b;
            state->particles[rand]->r = state->particles.pModel.r;
            state->particles[rand]->rf = state->particles.pModel.rf;

            //Set qDisp
            Eigen::Vector3d v = Random::get_vector();
            state->particles[rand]->qDisp = v;
            state->particles[rand]->qDisp = state->particles[rand]->qDisp.normalized() * state->particles[rand]->b;
            state->particles[rand]->pos = state->particles[rand]->com + state->particles[rand]->qDisp;

            state->particles[rand]->name = "Na";
            state->particles.cTot++;
            state->particles.aTot--;      
        }

        attempted++;
        state->move_callback({static_cast<unsigned int>(rand)});
    }

    bool accept(double dE) const{
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
        return Move::dump();
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
    mutable int pAtt = 0, nAtt = 0, pAcc = 0, nAcc = 0;
    mutable int* att;
    mutable int* acc;

    public:

    GrandCanonicalSingle(double chemPot, double donnan, double w) : Move(0.0, w),
                  d(donnan), cp(chemPot){
        if(ADD){
            this->id = "GCAddSingle";
        }
        else{
            this->id = "GCRemSingle";
        }
        
        printf("\t%s\n", this->id.c_str());
        constants::cp = chemPot;
        this->pVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.pModel.rf);
        this->nVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.nModel.rf);
        printf("\tCation accessible volume: %.3lf, Anion accessible volume: %.3lf\n", this->pVolume, this->nVolume);
        printf("\tChemical potential: %.3lf, Bias potential: %.3lf\n", this->cp, this->d);
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){
        if(ADD){
            auto [ind, qt] = state->particles.add_random(state->geo->_dh);
            this->q = qt;
            state->move_callback({ind});
        }

        else{
            std::vector< unsigned int > particles;
            auto [ind, qt] = state->particles.remove_random();
            this->q = qt;
            if(ind != -1){
                particles.push_back(ind);
            }
            state->move_callback(particles);     
        }
        this->attempted++;
    }

    bool accept(double dE) const{
        auto prob = -1.0;

        // ADD
        if(ADD){
            //Cation
            if(this->q > 0.0){
                prob = this->pVolume / state->particles.cTot * std::exp(this->cp - this->d * this->q - dE); //N + 1 since state->particles.cTot is the new N + 1 state
                this->acc = &this->pAcc;
                this->pAtt++;
            }
            //Anion
            else{
                prob = this->nVolume / state->particles.aTot * std::exp(this->cp - this->d * this->q - dE);
                this->acc = &this->nAcc;
                this->nAtt++;
            } 
        }
        // REMOVE
        else{
            //Cation
            if(this->q > 0.0){
                prob = (state->particles.cTot + 1) / this->pVolume * std::exp(this->d * this->q - this->cp - dE); //N since state->particles.cTot is the new N - 1 state
                this->acc = &this->pAcc;
                this->pAtt++;
            }
            //Anion
            else{
                prob = (state->particles.aTot + 1) / this->nVolume * std::exp(this->d * this->q - this->cp - dE);
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

    GrandCanonical(double chemPot, double w) : Move(0.0, w), cp(chemPot){
        if(ADD){
            this->id = "GCAdd";
        }
        else{
            this->id = "GCRem";
        }
        
        constants::cp = chemPot;
        this->pVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.pModel.rf);
        this->nVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.nModel.rf);
        this->valency = state->particles.pModel.q;

        printf("\t%s\n", this->id.c_str());
        printf("\tCation accessible volume: %.3lf, Anion accessible volume: %.3lf\n", this->pVolume, this->nVolume);
        printf("\tChemical potential: %.3lf", this->cp);
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){
        std::vector< unsigned int > particles;
        if(ADD){
            auto [ind, qt] = state->particles.add_random(state->geo->_dh, 1.0);
            particles.push_back(ind);
            for(int i = 0; i < this->valency; i++){
                auto [ind2, qt2] = state->particles.add_random(state->geo->_dh, -1.0);
                particles.push_back(ind2);
            }
        }
        else{
            std::vector< unsigned int > pt;

            int rand = 0;
            if(state->particles.cTot > 0){
                do{
                    rand = Random::get_random(state->particles.tot);
                } while(state->particles[rand]->q != state->particles.pModel.q);
                particles.push_back(rand);
            }

            if(state->particles.aTot > 0){
                while(particles.size() < this->valency + 1){
                    do{
                        rand = Random::get_random(state->particles.tot);
                    } while(state->particles[rand]->q != state->particles.nModel.q);

                    if(std::none_of(particles.begin(), particles.end(), [&](int val){return val==rand;})){
                        particles.push_back(rand);   
                    }   
                }
            }
            std::sort(particles.begin(), particles.end());
            std::reverse(particles.begin(), particles.end());
            for(auto p : particles){
                state->particles.remove(p);
            }
        }

        this->attempted++;
        state->move_callback(particles);  
    }

    bool accept(double dE) const{
        double prob;

        // ADD
        if(ADD){
            auto pFac = this->pVolume / state->particles.cTot;
            auto nFac = 1.0;
            for(int i = 0; i < this->valency; i++){
                nFac *= this->nVolume / (state->particles.aTot - i);
            }

            prob = pFac * nFac * std::exp((this->valency + 1.0) * this->cp - dE); 

        }
        // REMOVE
        else{
            auto pFac = (state->particles.cTot + 1) / this->pVolume;
            auto nFac = 1.0;
            for(int i = 0; i < this->valency; i++){
                nFac *= (state->particles.aTot + 1 - i) / this->nVolume ;
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
        return Move::dump();
    }
};


















class VolumeMove: public Move{
    private:
    double _oldV;
    double pressure;
    static constexpr auto unit = 2.430527863808942e-10;

    public:
    VolumeMove(double step, double pressure, double w) : Move(step, w){
        this->id = "Vol";
        printf("\t%s\n", this->id.c_str());
        this->pressure = pressure * unit;
        printf("\tPressure: %lf\n", this->pressure);
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        _oldV = state->geo->volume;
        auto lnV = std::log(state->geo->volume) + (Random::get_random() * 2.0 - 1.0) * this->stepSize;
        auto V = std::exp(lnV);
        auto L = std::cbrt(V);
        auto RL = L / state->geo->_d[0];

        std::vector<double> LV = {L, L, L};
        std::vector<double> LVh = {L / 2.0, L / 2.0, L / 2.0};

        state->geo->_d = LV;
        state->geo->_dh = LVh;
        state->geo->d = LV;
        state->geo->dh = LVh;
        state->geo->volume = V;

        std::vector< unsigned int > particles;

        for(unsigned int i = 0; i < state->particles.tot; i++){
            state->particles[i]->com *= RL;
            state->particles[i]->pos = state->particles[i]->com + state->particles[i]->qDisp;
            particles.push_back(state->particles[i]->index);
        }

        this->attempted++;
        state->move_callback(particles);
    }

    bool accept(double dE) const{
        auto ret = false;
        auto prob = exp(-dE - this->pressure * (state->geo->volume - _oldV) +
                      (state->particles.tot + 1.0) * std::log(state->geo->volume / _oldV));

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
        return Move::dump();
    }
};




class ChargeTrans: public Move{

    public:

    ChargeTrans(double step, double w) : Move(step, w){
        this->id = "qTrans";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        auto rand = 0;
        do{
            rand = Random::get_random(state->particles.tot);
        } while(state->particles[rand]->q < 0.0);
        
        state->particles[rand]->chargeTrans(this->stepSize);

        this->attempted++;
        state->move_callback({state->particles[rand]->index});
    }

    bool accept(double dE) const{
        auto ret = false;

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
        return Move::dump();
    }
};




class ChargeTranslate: public Move{

    public:

    ChargeTranslate(double step, double w) : Move(step, w){
        this->id = "qTranslate";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        auto rand = 0;
        do{
            rand = Random::get_random(state->particles.tot);
        } while(state->particles[rand]->q < 0.0);
        
        state->particles[rand]->chargeTranslate(this->stepSize);

        state->particles[rand]->qDisp = state->geo->displacement(state->particles[rand]->pos, state->particles[rand]->com);
        state->particles[rand]->b = state->particles[rand]->qDisp.norm();

        this->attempted++;
        state->move_callback({state->particles[rand]->index});
    }

    bool accept(double dE) const{
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
        return Move::dump();
    }
};



class ChargeTransRand: public Move{

    public:

    ChargeTransRand(double step, double w) : Move(step, w){
        this->id = "qTransRand";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        //printf("Move\n");
        std::vector< unsigned int > particles;
        auto rand = 0;
        if(state->particles.cTot > 0){
            do{
                rand = Random::get_random(state->particles.tot);
            } while(state->particles[rand]->q < 0.0);

            state->particles[rand]->chargeTransRand();

            particles.push_back(state->particles[rand]->index);
        }

        this->attempted++;
        state->move_callback(particles);
    }

    bool accept(double dE) const{
        auto ret = false;

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
        return Move::dump();
    }
};


class Cluster : public Move{
    private:
    double minDist;
    bool found = false;
    mutable std::unordered_map<unsigned int, unsigned int> att;
    mutable std::unordered_map<unsigned int, unsigned int> acc;
    unsigned int pNum = 0;
    std::shared_ptr<Particle> p;

    public:

    Cluster(double step, double md, double w) : Move(step, w), minDist(md){
        this->id = "Clus";
        printf("\t%s\n", this->id.c_str());
        printf("\tStepsize: %lf\n", step);
        printf("\tMax distance in cluster: %lf\n", this->minDist);
        printf("\tWeight: %lf\n", this->weight);
    }


    void operator()(){
        this->p = state->particles.random();
        std::vector<unsigned int> indices;
        Eigen::Vector3d disp;

        for(auto i : state->particles.particles){
            if(i->index == this->p->index) continue;

            if(state->geo->distance(i->pos, this->p->pos) <= this->minDist){
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
            state->particles.translate(indices, disp);
            state->move_callback(indices);
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

    bool accept(double dE) const{
        auto count = 1;

        if(found){

            for(auto i : state->particles.particles){
                if(i->index == this->p->index) continue;
                if(state->geo->distance(i->pos, this->p->pos) <= this->minDist){
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
    mutable double cp = 0.0;
    mutable int samples = 0;
    public:

    WidomInsertion(double w) : Move(0.0, w) {

        this->id = "WIns";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){      
        auto [ind, qt] = state->particles.add_random(state->geo->_dh, 2.0);
        this->q = qt;
        this->attempted++;
        state->move_callback({ind});
    }

    bool accept(double dE) const{
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
    mutable double cp = 0.0;
    mutable int samples = 0;
    public:

    WidomDeletion(double w) : Move(0.0, w) {

        this->id = "WDel";
        printf("\t%s\n", this->id.c_str());
        printf("\tWeight: %lf\n", this->weight);
    }

    void operator()(){      
        auto [ind, qt] = state->particles.remove_random();
        this->q = qt;

        this->attempted++;
        state->move_callback({ind});
    }

    bool accept(double dE) const{
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


Move* Move::create_move(MoveTypes moveType, std::vector<double> args){
    switch(moveType){
        case MoveTypes::Translate:
            assert(args.size() == 2);
                                        //step, w
            return new Translate(args[1], args[0]);
            break;
        case MoveTypes::GCSingleAdd:
            assert(args.size() == 3);
                                                            //cp, d, w
            return new GrandCanonicalSingle<true>(args[1], args[2], args[0]);
            break;
        case MoveTypes::GCSingleRemove:
            assert(args.size() == 3);
            return new GrandCanonicalSingle<false>(args[1], args[2], args[0]);
            break;
        case MoveTypes::Rotate:
            assert(args.size() == 2);
                                        //step, w
            return new Rotate(args[0], args[1]);
            break;
        case MoveTypes::Swap:
            assert(args.size() == 1);
            return new Swap(args[0]);
            break;
        case MoveTypes::SingleSwap:
            assert(args.size() == 1);
            return new SingleSwap(args[0]);
            break;
        case MoveTypes::VolumeMove:
            assert(args.size() == 3);
                                        //step, press, v, w
            return new VolumeMove(args[0], args[2], args[1]);
            break;
        case MoveTypes::ChargeTrans:
            assert(args.size() == 2);
            return new ChargeTrans(args[0], args[1]);
            break;
        case MoveTypes::ChargeTransRand:
            assert(args.size() == 2);
            return new ChargeTransRand(args[0], args[1]);
            break;
        case MoveTypes::Cluster:
            assert(args.size() == 3);
            return new Cluster(args[0], args[1], args[2]);
            break;
        case MoveTypes::WidomInsertion:
            assert(args.size() == 1);
            return new WidomInsertion(args[0]);
            break;
        case MoveTypes::WidomDeletion:
            assert(args.size() == 1);
            return new WidomDeletion(args[0]);
            break;
        case MoveTypes::GCAdd:
            assert(args.size() == 2);
                                                    //chempot, w
            return new GrandCanonical<true>(args[1], args[0]);
            break;
        case MoveTypes::GCRemove:
            assert(args.size() == 2);
            return new GrandCanonical<false>(args[1], args[0]);
            break;
        case MoveTypes::ChargeTranslate:
            assert(args.size() == 2);
            return new ChargeTranslate(args[0], args[1]);
            break;
        default:
            printf("Invalid move");
            break;
    }
}




}