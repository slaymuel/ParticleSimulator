#include "particles.h"

namespace Simulator{

bool Particles::is_valid(){
    return setPModel && setNModel && particles.size() > 0;
}

Eigen::MatrixXd Particles::get_subset(int sr, int fr) const{
    return this->positions.block(sr, 0, fr, 3);
}
Eigen::MatrixXd Particles::get_particle_pos(int index) const{
    return this->positions.row(index);
}

Particles::Particles(){}

std::shared_ptr<Particle> Particles::operator[](std::size_t index){
    return particles[index];
}

const std::shared_ptr<Particle> Particles::operator[](std::size_t index) const{
    return particles[index];
}


unsigned int Particles::size() const{
    return particles.size();
}

std::vector< std::shared_ptr<Particle> > Particles::get_subset(std::vector<unsigned int> &ps) const{
    // Should use range_view instead
    std::vector< std::shared_ptr<Particle> > subset;
    for (auto p : ps)
        if (p < this->tot) subset.push_back(this->particles[p]);

    return subset;
}


std::vector< std::shared_ptr<Particle> > Particles::get_subset(unsigned int i) const{
    std::vector< std::shared_ptr<Particle> > subset;
    subset.push_back(this->particles[i]);
    return subset;
}


std::shared_ptr<Particle> Particles::random() const{
    //std::sample
    return this->particles[Random::get_random(this->tot)];
}

void Particles::translate(std::vector<unsigned int> &ps, std::vector<double> &disp){
    for(auto p : ps){
        this->particles[p]->translate(disp);
    }
}

void Particles::translate(std::vector<unsigned int> &ps, Eigen::Vector3d &disp){
    for(auto p : ps){
        this->particles[p]->translate(disp);
    }
}


void Particles::add(const std::shared_ptr<Particle>& p, int index){
    int addAtIndex = index;
    if(index >= 0){
        this->particles.push_back(std::make_shared<Particle>());
        for(int i = this->tot; i > index; i--){
            *(this->particles[i]) = *(this->particles[i - 1]);
            this->particles[i]->index = i;
        }
    }
    else {
        addAtIndex = this->tot;
        this->particles.push_back(std::make_shared<Particle>());
    }
        
    this->particles[addAtIndex]->com << p->com[0], p->com[1], p->com[2];
    this->particles[addAtIndex]->pos << p->pos[0], p->pos[1], p->pos[2];
    this->particles[addAtIndex]->qDisp << p->qDisp;

    this->particles[addAtIndex]->index = addAtIndex;
    this->particles[addAtIndex]->r = p->r;
    this->particles[addAtIndex]->rf = p->rf;
    this->particles[addAtIndex]->q = p->q;
    this->particles[addAtIndex]->b = p->b;
    this->particles[addAtIndex]->b_min = p->b_min;
    this->particles[addAtIndex]->b_max = p->b_max;
    this->particles[addAtIndex]->name = p->name;

    if(p->q > 0){
        this->cTot++;
    }
    else {
        this->aTot++;
    }      
    this->tot++;
}



std::tuple<unsigned int, double> Particles::add_random(std::vector<double> box, double type){
    double rand = Random::get_random();
    double q;
    Eigen::Vector3d com;
    
    if(type != 0) rand = type;
    //Add cation
    if(rand >= 0.5){
        com = Random::random_pos_box(this->pModel.rf, box);
        //Eigen::Vector3d com;
        //com << 18.0, 23.0, 34.0;
        this->add(com, this->pModel.r, this->pModel.rf, this->pModel.q, this->pModel.b_min, this->pModel.b_max, "Na");
        q = this->pModel.q;
    }

    //Add anion
    else{
        com = Random::random_pos_box(this->nModel.rf, box);
        this->add(com, this->nModel.r, this->nModel.rf, this->nModel.q, this->nModel.b_min, this->nModel.b_max, "Cl");
        q = this->nModel.q;
    }
    //printf("Adding %i charge %lf com %lf %lf %lf\n", this->tot - 1, q, com[0], com[1], com[2]);
    return {this->tot - 1, q};
}



std::tuple<unsigned int, double> Particles::remove_random(){
    double q;
    double rand = Random::get_random();
    int rand2 = Random::get_random(this->tot);

    if(rand < 0.5){
        if(this->cTot > 0){
            do{
                rand2 = Random::get_random(this->tot);
            } while(this->particles[rand2]->q != this->pModel.q);
            q = this->particles[rand2]->q;
            this->remove(rand2);
        }

        else{
            rand2 = -1;
        }
    }
    else{
        if(this->aTot > 0){
            do{
                rand2 = Random::get_random(this->tot);
            } while(this->particles[rand2]->q != this->nModel.q);
            q = this->particles[rand2]->q;
            this->remove(rand2);
        }
        else{
            rand2 = -1;
        }
    }

    return {rand2, q};
}


std::tuple<unsigned int, double> Particles::remove_random(double type){
    double q;

    int rand = Random::get_random(this->tot);

    if(type > 0.0){
        if(this->cTot > 0){
            do{
                rand = Random::get_random(this->tot);
            } while(this->particles[rand]->q != this->pModel.q);
            q = this->particles[rand]->q;
            this->remove(rand);
        }
        else{
            rand = -1;
        }
    }
    else{
        if(this->aTot > 0){
            do{
                rand = Random::get_random(this->tot);
            } while(this->particles[rand]->q != this->nModel.q);
            q = this->particles[rand]->q;
            this->remove(rand);
        }
        else{
            rand = -1;
        }
    }
    return {rand, q};
}


void Particles::remove(std::size_t index){
    if(this->particles[index]->q > 0){
        this->cTot--;
    }
    else{
        this->aTot--; 
    }
    
    //move last particle in particles to position of particle to be removed
    this->particles.erase(this->particles.begin() + index);

    for(unsigned int i = index; i < this->tot - 1; i++)
        this->particles[i]->index = i;

    this->tot--;
}


void Particles::create(int pNum, int nNum, std::map<std::string, double> params){

    this->pModel.q = params["p"];
    this->pModel.r = params["rp"];
    this->pModel.rf = params["rfp"];
    this->pModel.b_min = params["bp_min"];
    this->pModel.b_max = params["bp_max"];
    this->pModel.name = "Na";

    this->nModel.q = params["n"];
    this->nModel.r = params["rn"];
    this->nModel.rf = params["rfn"];
    this->nModel.b_min = params["bn_min"];
    this->nModel.b_max = params["bn_max"];
    this->nModel.name = "Cl";

    Eigen::Vector3d com;
    for(int i = 0; i < pNum + nNum; i++){
        com = Random::get_vector();
        (i < pNum) ? this->add(com, this->pModel.r, this->pModel.rf, this->pModel.q, this->pModel.b_min, this->pModel.b_max, "Na") : 
                        this->add(com, this->nModel.r, this->nModel.rf, this->nModel.q, this->nModel.b_min, this->nModel.b_max, "Cl");
    }

    printf("\n");
    Logger::Log("Created ", pNum, " cations and ", nNum, " anions");
    for(auto val : params){
        Logger::Log("\t", val.first.c_str(), " ", val.second);
    }
    printf("\n");
    if(params.size() != 10){
        printf("Wrong number of arguments in particles::create params, contains %lu\n", params.size());
        exit(1);
    }
}


} // end of namespace Simulator