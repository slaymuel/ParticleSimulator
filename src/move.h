#pragma once
#include "particle.h"
#include "particles.h"
#include <math.h>
#include "state.h"
#include <algorithm>
#include <unordered_map>
#include "logger.h"
#include "genericfactory.h"

namespace Simulator{

namespace Moves{


class Move;

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

// Move factory
using moveCreator = std::function<std::unique_ptr<Move>(std::string, std::vector<double>)>;
inline GenericFactory< Move, MoveTypes, moveCreator> moveFactory;

// Base class for moves
class Move{

    protected:
    mutable int accepted, rejected, attempted;
    double stepSize;
    std::string id;

    public:
    double weight;
    static State* state;

    Move(double step, double w, std::string id);
    virtual ~Move() = default;

    virtual void operator()() = 0;
    virtual bool accept(double dE) const = 0;
    virtual std::string dump() const = 0;
};

//std::string Move::dump() const;

class Translate : public Move{
    public:

    Translate(double step, double w);
    ~Translate() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};



class Rotate : public Move{
    public:

    Rotate(double step, double w);
    ~Rotate() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};



class Swap : public Move{

    public:
    Swap(double w);
    ~Swap() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};





class SingleSwap : public Move{

    public:
    SingleSwap(double w);
    ~SingleSwap() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
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

namespace {
    std::unique_ptr<Move> createGrandCanonicalSingleADD(std::string name, std::vector<double> args){
        return std::make_unique< GrandCanonicalSingle<true> >(args[1], args[2], args[0]);
    }

    bool registered_grandCanonicalSingle_add = 
                    moveFactory.registerObject(MoveTypes::GCSINGLEADD, createGrandCanonicalSingleADD);
}
namespace {
    std::unique_ptr<Move> createGrandCanonicalSingleRemove(std::string name, std::vector<double> args){
        return std::make_unique< GrandCanonicalSingle<false> >(args[1], args[2], args[0]);
    }

    bool registered_grandCanonicalSingle_rem = 
                    moveFactory.registerObject(MoveTypes::GCSINGLEREMOVE, createGrandCanonicalSingleRemove);
}





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


namespace {
    std::unique_ptr<Move> createGrandCanonicalAdd(std::string name, std::vector<double> args){
        return std::make_unique< GrandCanonical<true> >(args[1], args[0]);
    }

    bool registered_GrandCanonicalADD = moveFactory.registerObject(MoveTypes::GCADD, createGrandCanonicalAdd);
}
namespace {
    std::unique_ptr<Move> createGrandCanonicalRem(std::string name, std::vector<double> args){
        return std::make_unique< GrandCanonical<false> >(args[1], args[0]);
    }

    bool registered_GrandCanonicalRem = moveFactory.registerObject(MoveTypes::GCREMOVE, createGrandCanonicalRem);
}










class VolumeMove: public Move{
    private:
    double _oldV;
    double pressure;
    static constexpr auto unit = 2.430527863808942e-10;

    public:
    VolumeMove(double step, double pressure, double w);
    ~VolumeMove() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};




class ChargeTrans: public Move{

    public:

    ChargeTrans(double step, double w);
    ~ChargeTrans() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};




class ChargeTranslate: public Move{

    public:

    ChargeTranslate(double step, double w);
    ~ChargeTranslate() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};



class ChargeTransRand: public Move{

    public:

    ChargeTransRand(double step, double w);
    ~ChargeTransRand() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
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

    Cluster(double step, double md, double w);
    ~Cluster() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};


class WidomInsertion : public Move{

    protected:
    double q;
    mutable double cp = 0.0;
    mutable int samples = 0;
    public:

    WidomInsertion(double w);
    ~WidomInsertion() override = default;

    void operator()() override;

    bool accept(double dE) const override;

    std::string dump() const override;
};

class WidomDeletion : public Move{

    protected:
    double q;
    mutable double cp = 0.0;
    mutable int samples = 0;
    public:

    WidomDeletion(double w);
    ~WidomDeletion() override = default;

    void operator()() override;
    bool accept(double dE) const override;
    std::string dump() const override;
};

} // end of namespace Moves


} // end of namespace Simulator