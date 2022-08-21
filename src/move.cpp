#include "move.h"

namespace Simulator{

namespace Moves{


Move::Move(double step, double w, std::string id) : stepSize(step), id(id), weight(w){
    Logger::Log("\t", this->id);
    Logger::Log("\tWeight: ", this->weight);
    this->accepted = 0;
    this->rejected = 0;
    this->attempted = 0;
}

std::string Move::dump() const{
    std::ostringstream s;
    s.precision(1);
    s << std::fixed;
    s << "\t" << this->id << ": " << (double) this->accepted / this->attempted * 100.0 << "%, " << 
                                                this->attempted << " (" << this->accepted <<") ";
    return s.str();
}

State* Move::state = nullptr;


Translate::Translate(double step, double w) : Move(step, w, "TRANSLATION"){
    Logger::Log("\tStepsize: ", step);
}

void Translate::operator()(){
    auto p = state->particles.random();

    p->translate(stepSize);
    this->attempted++;

    state->move_callback({p->index});
}

bool Translate::accept(double dE) const {
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
std::string Translate::dump() const {
    return Move::dump();
}

// Register the move to the move factory
namespace {
    std::unique_ptr<Move> createTranslate(std::string name, std::vector<double> args){
        return std::make_unique<Translate>(args[1], args[0]);
    }

    bool registered_translate = moveFactory.registerObject(MoveTypes::TRANSLATE, createTranslate);
}

       
            
            



Rotate::Rotate(double step, double w) : Move(step, w, "ROTATION"){
    Logger::Log("\tStepsize: ", step);
}

void Rotate::operator()() {
    auto p = state->particles.random();

    p->rotate(this->stepSize);
    attempted++;

    state->move_callback({p->index});
}

bool Rotate::accept(double dE) const {
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

std::string Rotate::dump() const {
    return Move::dump();
}

namespace {
    std::unique_ptr<Move> createRotate(std::string name, std::vector<double> args){
        return std::make_unique<Rotate>(args[0], args[1]);
    }

    bool registered_rotate = moveFactory.registerObject(MoveTypes::ROTATE, createRotate);
}


Swap::Swap(double w) : Move(0.0, w, "SWAP"){}

void Swap::operator()()  {
    int rand2 = -1;
    auto rand = Random::get_random(state->particles.tot);

    do{
        rand2 = Random::get_random(state->particles.tot);
    } while(state->particles[rand]->q == state->particles[rand2]->q);

    std::swap(state->particles[rand]->pos, state->particles[rand2]->pos);
    std::swap(state->particles[rand]->com, state->particles[rand2]->com);

    attempted++;
    state->move_callback({static_cast<unsigned int>(rand), static_cast<unsigned int>(rand2)});
}

bool Swap::accept(double dE) const {
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

std::string Swap::dump() const {
    return Move::dump();
}

namespace {
    std::unique_ptr<Move> createSwap(std::string name, std::vector<double> args){
        return std::make_unique<Swap>(args[0]);
    }

    bool registered_Swap = moveFactory.registerObject(MoveTypes::SWAP, createSwap);
}






SingleSwap::SingleSwap(double w) : Move(0.0, w, "SINGLESWAP"){}

void SingleSwap::operator()()  {
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

bool SingleSwap::accept(double dE) const {
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

std::string SingleSwap::dump() const  {
    return Move::dump();
}

namespace {
    std::unique_ptr<Move> createSingleSwap(std::string name, std::vector<double> args){
        return std::make_unique<SingleSwap>(args[0]);
    }

    bool registered_singleSwap = moveFactory.registerObject(MoveTypes::SINGLESWAP, createSingleSwap);
}




VolumeMove::VolumeMove(double step, double pressure, double w) : Move(step, w, "VOLUMEMOVE"){
    this->pressure = pressure * unit;
    Logger::Log("\tPressure: ", this->pressure);
    Logger::Log("\tStepsize: ", step);
}

void VolumeMove::operator()()  {
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

bool VolumeMove::accept(double dE) const {
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

std::string VolumeMove::dump() const {
    return Move::dump();
}

namespace {
    std::unique_ptr<Move> createVolumeMove(std::string name, std::vector<double> args){
        return std::make_unique<VolumeMove>(args[0], args[2], args[1]);
    }

    bool registered_VolumeMove = moveFactory.registerObject(MoveTypes::VOLUMEMOVE, createVolumeMove);
}


ChargeTrans::ChargeTrans(double step, double w) : Move(step, w, "CHARGETRANS"){
    Logger::Log("\tStepsize: ", step);
}

void ChargeTrans::operator()()  {
    auto rand = 0;
    do{
        rand = Random::get_random(state->particles.tot);
    } while(state->particles[rand]->q < 0.0);
    
    state->particles[rand]->chargeTrans(this->stepSize);

    this->attempted++;
    state->move_callback({state->particles[rand]->index});
}

bool ChargeTrans::accept(double dE) const {
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

std::string ChargeTrans::dump() const  {
    return Move::dump();
}

namespace {
    std::unique_ptr<Move> createChargeTrans(std::string name, std::vector<double> args){
        return std::make_unique<ChargeTrans>(args[0], args[1]);
    }

    bool registered_ChargeTrans = moveFactory.registerObject(MoveTypes::CHARGETRANS, createChargeTrans);
}


ChargeTranslate::ChargeTranslate(double step, double w) : Move(step, w, "CHARGETRANS"){
    Logger::Log("\tStepsize: ", step);
}


void ChargeTranslate::operator()()  {
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

bool ChargeTranslate::accept(double dE) const {
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

std::string ChargeTranslate::dump() const {
    return Move::dump();
}



ChargeTransRand::ChargeTransRand(double step, double w) : Move(step, w, "CHARGETRANSRAND"){}

void ChargeTransRand::operator()() {
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

bool ChargeTransRand::accept(double dE) const {
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

std::string ChargeTransRand::dump() const {
    return Move::dump();
}

namespace {
    std::unique_ptr<Move> createChargeTransRand(std::string name, std::vector<double> args){
        return std::make_unique<ChargeTransRand>(args[0], args[1]);
    }

    bool registered_ChargeTransRand = moveFactory.registerObject(MoveTypes::CHARGETRANSRAND, createChargeTransRand);
}


Cluster::Cluster(double step, double md, double w) : Move(step, w, "CLUSTER"), minDist(md){
    Logger::Log("\tStepsize: ", step);
    Logger::Log("\tMax distance in cluster: ", this->minDist);
}

void Cluster::operator()() {
    this->p = state->particles.random();
    std::vector<unsigned int> indices;
    Eigen::Vector3d disp;

    for(auto i : state->particles){
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

bool Cluster::accept(double dE) const {
    auto count = 1;

    if(found){

        for(auto i : state->particles){
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

std::string Cluster::dump() const {
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

namespace {
    std::unique_ptr<Move> createCluster(std::string name, std::vector<double> args){
        return std::make_unique<Cluster>(args[0], args[1], args[2]);
    }

    bool registered_Cluster = moveFactory.registerObject(MoveTypes::CLUSTER, createCluster);
}

WidomInsertion::WidomInsertion(double w) : Move(0.0, w, "WIDOMINSERTION") {}

void WidomInsertion::operator()() {      
    auto [ind, qt] = state->particles.add_random(state->geo->_dh, 2.0);
    this->q = qt;
    this->attempted++;
    state->move_callback({ind});
}

bool WidomInsertion::accept(double dE) const {
    if(this->samples > 10000){
        this->samples = 0;
        this->cp = 0.0;
    }
    
    this->cp += std::exp(-dE);
    this->samples++;

    return false;
}

std::string WidomInsertion::dump() const {
    std::ostringstream ss;
    ss.precision(5);
    ss << std::fixed;
    ss << "\t" << this->id << " cp: " << -std::log(this->cp / this->samples) <<" " << this->cp / this->samples << " " << this->cp << " samples: " << this->samples << " attempted: " << this->attempted;
    return ss.str();
}


namespace {
    std::unique_ptr<Move> createWidomInsertion(std::string name, std::vector<double> args){
        return std::make_unique<WidomInsertion>(args[0]);
    }

    bool registered_WidomInsertion = moveFactory.registerObject(MoveTypes::WIDOMINSERTION, createWidomInsertion);
}


WidomDeletion::WidomDeletion(double w) : Move(0.0, w, "WIDOMDELETION") {}

void WidomDeletion::operator()() {      
    auto [ind, qt] = state->particles.remove_random();
    this->q = qt;

    this->attempted++;
    state->move_callback({ind});
}

bool WidomDeletion::accept(double dE) const {
    if(this->samples > 10000){
        this->samples = 0;
        this->cp = 0.0;
    }
    this->cp += std::exp(dE);
    this->samples++;

    return false;
}

std::string WidomDeletion::dump() const {
    std::ostringstream ss;
    ss.precision(5);
    ss << std::fixed;
    ss << "\t" << this->id << " cp: " << -std::log(this->cp / this->samples) <<" " << this->cp / this->samples << " " << this->cp << " samples: " << this->samples << " attempted: " << this->attempted;
    return ss.str();
}

namespace {
    std::unique_ptr<Move> createWidomDeletion(std::string name, std::vector<double> args){
        return std::make_unique<WidomDeletion>(args[0]);
    }

    bool registered_WidomDeletion = moveFactory.registerObject(MoveTypes::WIDOMDELETION, createWidomDeletion);
}


} // end of namespace Moves


} // end of namespace Simulator