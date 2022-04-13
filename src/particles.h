#pragma once

#include <tuple>
#include <map>

#include "particle.h"

namespace Simulator{

class Particles{
    private:

    public:
    bool setPModel = false;
    bool setNModel = false;
    Eigen::MatrixXd positions;
    Particle pModel;
    Particle nModel;
    std::vector< std::shared_ptr<Particle> > particles, cations, anions;
    std::vector<int> movedParticles;
    unsigned int cTot = 0, aTot = 0, tot = 0;

    Eigen::MatrixXd get_subset(int sr, int fr) const{
        return this->positions.block(sr, 0, fr, 3);
    }
    Eigen::MatrixXd get_particle_pos(int index) const{
        return this->positions.row(index);
    }

    Particles(){}

    std::shared_ptr<Particle> operator[](std::size_t index){
        return particles[index];
    }

    const std::shared_ptr<Particle> operator[](std::size_t index) const{
        return particles[index];
    }


    unsigned int size() const{
        return particles.size();
    }

    std::vector< std::shared_ptr<Particle> > get_subset(std::vector<unsigned int> &ps) const{
        std::vector< std::shared_ptr<Particle> > subset;

        //std::vector<Particle> subset(ps.size(), 0);
        //std::transform(ps.begin(), ps.end(), subset.begin(), [particles](size_t i) {return particles[i];});
        for (auto p : ps){
            if (p < this->tot) subset.push_back(this->particles[p]);
        }
        return subset;
    }


    std::vector< std::shared_ptr<Particle> > get_subset(unsigned int i) const{
        std::vector< std::shared_ptr<Particle> > subset;
        subset.push_back(this->particles[i]);
        return subset;
    }


    std::shared_ptr<Particle> random() const{
        //https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
        //std::sample

        //return this->particles[(*distribution)(rand_gen)];
        return this->particles[Random::get_random(this->tot)];
    }

    void translate(std::vector<unsigned int> &ps, std::vector<double> &disp){
        for(auto p : ps){
            this->particles[p]->translate(disp);
        }
    }

    void translate(std::vector<unsigned int> &ps, Eigen::Vector3d &disp){
        for(auto p : ps){
            this->particles[p]->translate(disp);
        }
    }




    template <typename T, typename G>
    void add(T com, T pos, G qDisp, double r, double rf, double q, double b, double b_min, double b_max, std::string name){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[this->tot]->com << com[0], com[1], com[2];
        this->particles[this->tot]->pos << pos[0], pos[1], pos[2];
        this->particles[this->tot]->qDisp << qDisp[0], qDisp[1], qDisp[2];

        if(std::abs(this->particles[this->tot]->qDisp.norm() - b) > 1e-10){
            printf("Error reading particles, b is not equal to |qDisp| for particle %u\n", this->tot);
            printf("b = %.10lf\n", b);
            printf("|qDisp| = %.10lf\n", this->particles[this->tot]->qDisp.norm());
            exit(1);   
        }

        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->b_min = b_min;
        this->particles[this->tot]->b_max = b_max;
        this->particles[this->tot]->name = name;
            
        if(q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->tot++;
        assert(this->tot <= this->particles.size() && "tot is larger than particle vector size\n");

        if(!this->setPModel){
            if(q > 0){
                this->pModel.q = q;
                this->pModel.b_min = b_min;
                this->pModel.b_max = b_max;
                this->pModel.r = r;
                this->pModel.rf = rf;
                this->pModel.name = name;
                this->setPModel = true;
            }
        }

        if(!this->setNModel){
            if(q < 0){
                this->nModel.q = q;
                this->nModel.b_min = b_min;
                this->nModel.b_max = b_max;
                this->nModel.r = r;
                this->nModel.rf = rf;
                this->nModel.name = name;
                this->setNModel = true;
            }
        }
    }



    template <typename T>
    void add(T com, double r, double rf, double q, double b_min, double b_max, std::string name){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);

        this->particles[this->tot]->com = com;
        
        
        //this->particles[this->tot]->qDisp = Random::get_norm_vector();

        //this->particles[this->tot]->b = b_min + (b_max - b_min) * Random::get_random();
        //this->particles[this->tot]->qDisp = Random::get_random_vector(1.0);
        //this->particles[this->tot]->qDisp = this->particles[this->tot]->qDisp.normalized() * this->particles[this->tot]->b;

        this->particles[this->tot]->qDisp = Random::get_random_vector(b_max);
        this->particles[this->tot]->b = this->particles[this->tot]->qDisp.norm();
        //std::cout << this->particles[this->tot]->b << std::endl;
        /*if(q > 0){
            this->particles[this->tot]->b = Random::get_random_from_distribution();
        }
        else{
            this->particles[this->tot]->b = 0.0;
        }*/

        this->particles[this->tot]->pos = this->particles[this->tot]->com + this->particles[this->tot]->qDisp;
        

        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b_min = b_min;
        this->particles[this->tot]->b_max = b_max;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->name = name;
        this->particles[this->tot]->index = this->tot;


        if(q > 0){
            this->cTot++;
        }
        else{
            this->aTot++;
        }

        this->tot++;
        assert(this->tot <= this->particles.size() && "tot is larger than particle vector size\n");
    }


    void add(const std::shared_ptr<Particle>& p){

        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}

        this->particles[this->tot]->com << p->com[0], p->com[1], p->com[2];
        this->particles[this->tot]->pos << p->pos[0], p->pos[1], p->pos[2];
        this->particles[this->tot]->qDisp << p->qDisp;

        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = p->r;
        this->particles[this->tot]->rf = p->rf;
        this->particles[this->tot]->q = p->q;
        this->particles[this->tot]->b = p->b;
        this->particles[this->tot]->b_min = p->b_min;
        this->particles[this->tot]->b_max = p->b_max;
        this->particles[this->tot]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->tot++;
    }


    void add(const std::shared_ptr<Particle>& p, int index){

        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}
        //std::copy(this->particles.begin() + index, this->particles.begin() + this->tot, this->particles.begin() + index + 1);
        for(int i = this->tot; i > index; i--){
            //printf("Moving particle %i to %i\n", i - 1, i);
            *(this->particles[i]) = *(this->particles[i - 1]);
            this->particles[i]->index = i;
        }

        this->tot++;

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[index]->com << p->com[0], p->com[1], p->com[2];
        this->particles[index]->pos << p->pos[0], p->pos[1], p->pos[2];
        this->particles[index]->qDisp << p->qDisp;

        this->particles[index]->index = index;
        this->particles[index]->r = p->r;
        this->particles[index]->rf = p->rf;
        this->particles[index]->q = p->q;
        this->particles[index]->b = p->b;
        this->particles[index]->b_min = p->b_min;
        this->particles[index]->b_max = p->b_max;
        this->particles[index]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }      
    }



    std::tuple<unsigned int, double> add_random(std::vector<double> box, double type = 0.0){
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



    std::tuple<unsigned int, double> remove_random(){
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

        /*if(this->cTot <= 0 || this->aTot <= 0){
            printf("Trying to remove last particle!\n");
        }*/

        //printf("Removing %i charge %lf\n", rand2, q);
        return {rand2, q};
    }


    std::tuple<unsigned int, double> remove_random(double type){
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


    void remove(std::size_t index){
        //this->particles.erase(this->particles.begin() + i);
        //printf("Copying in remove\n");
        //std::copy(this->particles.begin() + i + 1, this->particles.begin() + this->tot, this->particles.begin() + i);
        //std::copy(this->particles.begin() + i, this->particles.begin() + this->tot, this->particles.begin() + i - 1);

        if(this->particles[index]->q > 0){
            this->cTot--;
        }
        else{
            this->aTot--; 
        }
        
        //move last particle in particles to position of particle to be removed
        this->particles.erase(this->particles.begin() + index);

        for(unsigned int i = index; i < this->tot - 1; i++){
            //printf("Remove: Moving particle %i to %i\n", i + 1, i);
            //*(this->particles[i]) = *(this->particles[i + 1]);
            this->particles[i]->index = i;
        }

        this->tot--;
    }


    //void create(int pNum, int nNum, double p, double n, double rfp = 2.5, double rfn = 2.5, double rp = 2.5, double rn = 2.5, double bp = 0.0, double bp_min = 0.0, double bn = 0.0, double bn_min = 0.0){
    void create(int pNum, int nNum, std::map<std::string, double> params){

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
        printf("\n");
        for(auto val : params){
            Logger::Log("\t", val.first.c_str(), " ", val.second);
            printf("\n");
        }

        if(params.size() != 10){
            printf("Wrong number of arguments in particles::create params, contains %lu\n", params.size());
            exit(1);
        }
    }
};

}