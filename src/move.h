#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"
#include <algorithm>
#include <unordered_map>
#include "logger.h"

namespace Simulator{

using CallBack = std::function<void(std::vector< unsigned int >)>;

enum class MoveTypes{
    TRANSLATE,
    GCSINGLEADD,
    GCSINGLEREMOVE,
    ROTATE,
    SWAP,
    SINGLESWAP,
    VOLUMEMOVE,
    CHARGETRANS,
    CHARGETRANSRAND,
    CLUSTER,
    WIDOMINSERTION,
    WIDOMDELETION,
    GCADD,
    GCREMOVE,
    CHARGETRANSLATE
};




class Move{

    protected:
    mutable int accepted, rejected, attempted;
    double stepSize;
    std::string id;

    public:
    double weight;
    static State* state;

    Move(double step, double w, std::string id) : stepSize(step), id(id), weight(w){
        Logger::Log("\t", this->id);
        Logger::Log("\tWeight: ", this->weight);
        this->accepted = 0;
        this->rejected = 0;
        this->attempted = 0;
    }
    virtual ~Move() = default;

    virtual void operator()() = 0;
    virtual bool accept(double dE) const = 0;
    virtual std::string dump() const = 0;

    static std::unique_ptr<Move> createMove(MoveTypes move_type, std::vector<double> args);
};

std::string Move::dump() const{
    std::ostringstream s;
    s.precision(1);
    s << std::fixed;
    s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << this->attempted << " (" << this->accepted <<") ";
    return s.str();
}

State* Move::state = nullptr;


class Translate : public Move{
    public:

    Translate(double step, double w) : Move(step, w, "TRANSLATION"){
        Logger::Log("\tStepsize: ", step);
    }
    ~Translate() override = default;

    void operator()() override{
        auto p = state->particles.random();

        p->translate(stepSize);
        this->attempted++;

        state->move_callback({p->index});
    }

    bool accept(double dE) const override{
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
    std::string dump() const override{
        return Move::dump();
    }
};



class Rotate : public Move{
    public:

    Rotate(double step, double w) : Move(step, w, "ROTATION"){
        Logger::Log("\tStepsize: ", step);
    }
    ~Rotate() override = default;

    void operator()() override{
        auto p = state->particles.random();

        p->rotate(this->stepSize);
        attempted++;

        state->move_callback({p->index});
    }

    bool accept(double dE) const override{
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
    std::string dump() const override{
        return Move::dump();
    }
};



class Swap : public Move{

    public:
    Swap(double w) : Move(0.0, w, "SWAP"){}
    ~Swap() override = default;

    void operator()() override {
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

    bool accept(double dE) const override{
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

    std::string dump() const override{
        return Move::dump();
    }
};





class SingleSwap : public Move{

    public:
    SingleSwap(double w) : Move(0.0, w, "SINGLESWAP"){}
    ~SingleSwap() override = default;

    void operator()() override {
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

    bool accept(double dE) const override{
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

    std::string dump() const override {
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

    GrandCanonicalSingle(double chemPot, double donnan, double w) : Move(0.0, w, ADD ? "GCSINGLEADD" : "GCSINGLEREM"),
                  d(donnan), cp(chemPot){

        constants::cp = chemPot;
        this->pVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.pModel.rf);
        this->nVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.nModel.rf);
        Logger::Log("\tCation accessible volume: ", this->pVolume, " Anion accessible volume: ", this->nVolume);
        Logger::Log("\tChemical potential: ", this->cp, " Bias potential: ", this->d);
    }
    ~GrandCanonicalSingle() override = default;

    void operator()() override {
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

    bool accept(double dE) const override{
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

    std::string dump() const override{
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

    GrandCanonical(double chemPot, double w) : Move(0.0, w, ADD ? "GCADD" : "GCREM"), cp(chemPot){
        constants::cp = chemPot;
        this->pVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.pModel.rf);
        this->nVolume = state->geo->_d[0] * state->geo->_d[1] * (state->geo->_d[2] - 2.0 * state->particles.nModel.rf);
        this->valency = state->particles.pModel.q;

        Logger::Log("\tCation accessible volume: ", this->pVolume, ", Anion accessible volume: ", this->nVolume);
        Logger::Log("\tChemical potential: ", this->cp);
    }
    ~GrandCanonical() override = default;

    void operator()() override {
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

    bool accept(double dE) const override{
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

    std::string dump() const override {
        return Move::dump();
    }
};


















class VolumeMove: public Move{
    private:
    double _oldV;
    double pressure;
    static constexpr auto unit = 2.430527863808942e-10;

    public:
    VolumeMove(double step, double pressure, double w) : Move(step, w, "VOLUMEMOVE"){
        this->pressure = pressure * unit;
        Logger::Log("\tPressure: ", this->pressure);
        Logger::Log("\tStepsize: ", step);
    }
    ~VolumeMove() override = default;

    void operator()() override {
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

    bool accept(double dE) const override{
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

    std::string dump() const override{
        return Move::dump();
    }
};




class ChargeTrans: public Move{

    public:

    ChargeTrans(double step, double w) : Move(step, w, "CHARGETRANS"){
        Logger::Log("\tStepsize: ", step);
    }
    ~ChargeTrans() override = default;

    void operator()() override {
        auto rand = 0;
        do{
            rand = Random::get_random(state->particles.tot);
        } while(state->particles[rand]->q < 0.0);
        
        state->particles[rand]->chargeTrans(this->stepSize);

        this->attempted++;
        state->move_callback({state->particles[rand]->index});
    }

    bool accept(double dE) const override{
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

    std::string dump() const override {
        return Move::dump();
    }
};




class ChargeTranslate: public Move{

    public:

    ChargeTranslate(double step, double w) : Move(step, w, "CHARGETRANS"){
        Logger::Log("\tStepsize: ", step);
    }
    ~ChargeTranslate() override = default;

    void operator()() override {
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

    bool accept(double dE) const override{
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

    std::string dump() const override{
        return Move::dump();
    }
};



class ChargeTransRand: public Move{

    public:

    ChargeTransRand(double step, double w) : Move(step, w, "CHARGETRANSRAND"){}
    ~ChargeTransRand() override = default;

    void operator()() override{
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

    bool accept(double dE) const override{
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

    std::string dump() const override{
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

    Cluster(double step, double md, double w) : Move(step, w, "CLUSTER"), minDist(md){
        Logger::Log("\tStepsize: ", step);
        Logger::Log("\tMax distance in cluster: ", this->minDist);
    }
    ~Cluster() override = default;

    void operator()() override{
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

    bool accept(double dE) const override{
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

    std::string dump() const override{
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

    WidomInsertion(double w) : Move(0.0, w, "WIDOMINSERTION") {}
    ~WidomInsertion() override = default;

    void operator()() override{      
        auto [ind, qt] = state->particles.add_random(state->geo->_dh, 2.0);
        this->q = qt;
        this->attempted++;
        state->move_callback({ind});
    }

    bool accept(double dE) const override{
        if(this->samples > 10000){
            this->samples = 0;
            this->cp = 0.0;
        }
        
        this->cp += std::exp(-dE);
        this->samples++;

        return false;
    }

    std::string dump() const override{
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

    WidomDeletion(double w) : Move(0.0, w, "WIDOMDELETION") {}
    ~WidomDeletion() override = default;

    void operator()() override{      
        auto [ind, qt] = state->particles.remove_random();
        this->q = qt;

        this->attempted++;
        state->move_callback({ind});
    }

    bool accept(double dE) const override{
        if(this->samples > 10000){
            this->samples = 0;
            this->cp = 0.0;
        }
        this->cp += std::exp(dE);
        this->samples++;

        return false;
    }

    std::string dump() const override{
        std::ostringstream ss;
        ss.precision(5);
        ss << std::fixed;
        ss << "\t" << this->id << " cp: " << -std::log(this->cp / this->samples) <<" " << this->cp / this->samples << " " << this->cp << " samples: " << this->samples << " attempted: " << this->attempted;
        return ss.str();
    }
};


std::unique_ptr<Move> Move::createMove(MoveTypes moveType, std::vector<double> args){
    printf("\n");
    Logger::Log("Adding move:");
    switch(moveType){
        case MoveTypes::TRANSLATE:
            assert(args.size() == 2);
                                        //step, w
            return std::make_unique<Translate>(args[1], args[0]);
            break;
        case MoveTypes::GCSINGLEADD:
            assert(args.size() == 3);
                                                            //cp, d, w
            return std::make_unique< GrandCanonicalSingle<true> >(args[1], args[2], args[0]);
            break;
        case MoveTypes::GCSINGLEREMOVE:
            assert(args.size() == 3);
            return std::make_unique< GrandCanonicalSingle<false> >(args[1], args[2], args[0]);
            break;
        case MoveTypes::ROTATE:
            assert(args.size() == 2);
                                        //step, w
            return std::make_unique<Rotate>(args[0], args[1]);
            break;
        case MoveTypes::SWAP:
            assert(args.size() == 1);
            return std::make_unique<Swap>(args[0]);
            break;
        case MoveTypes::SINGLESWAP:
            assert(args.size() == 1);
            return std::make_unique<SingleSwap>(args[0]);
            break;
        case MoveTypes::VOLUMEMOVE:
            assert(args.size() == 3);
                                        //step, press, v, w
            return std::make_unique<VolumeMove>(args[0], args[2], args[1]);
            break;
        case MoveTypes::CHARGETRANS:
            assert(args.size() == 2);
            return std::make_unique<ChargeTrans>(args[0], args[1]);
            break;
        case MoveTypes::CHARGETRANSRAND:
            assert(args.size() == 2);
            return std::make_unique<ChargeTransRand>(args[0], args[1]);
            break;
        case MoveTypes::CLUSTER:
            assert(args.size() == 3);
            return std::make_unique<Cluster>(args[0], args[1], args[2]);
            break;
        case MoveTypes::WIDOMINSERTION:
            assert(args.size() == 1);
            return std::make_unique<WidomInsertion>(args[0]);
            break;
        case MoveTypes::WIDOMDELETION:
            assert(args.size() == 1);
            return std::make_unique<WidomDeletion>(args[0]);
            break;
        case MoveTypes::GCADD:
            assert(args.size() == 2);
                                                    //chempot, w
            return std::make_unique< GrandCanonical<true> >(args[1], args[0]);
            break;
        case MoveTypes::GCREMOVE:
            assert(args.size() == 2);
            return std::make_unique< GrandCanonical<false> >(args[1], args[0]);
            break;
        case MoveTypes::CHARGETRANSLATE:
            assert(args.size() == 2);
            return std::make_unique<ChargeTranslate>(args[0], args[1]);
            break;
        default:
            Logger::Log("Invalid move");
            break;
    }
}




}