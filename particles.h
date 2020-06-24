#pragma once

#include <Eigen/Dense>
#include <vector>
#include "particle.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include "../libxdrfile/include/xdrfile_xtc.h"

class Particles{

    public:
    //Eigen::MatrixXd positions;
    Particle pModel;
    Particle nModel;
    std::vector< std::shared_ptr<Particle> > particles, cations, anions;
    std::vector<int> movedParticles;
    unsigned int pTot = 0, cTot = 0, aTot = 0, iTot = 0, tot = 0;

    //Eigen::MatrixXd get_subset(int sr, int fr){
    //    return this->positions.block(sr, 0, fr, 3);
    //}

    Particles(){}

    std::shared_ptr<Particle> operator[](std::size_t index){
        return particles[index];
    }

    std::vector<std::shared_ptr<Particle>> get_subset(std::vector<unsigned int> &ps){
        std::vector<std::shared_ptr<Particle>> subset;
        for (auto p : ps){
            if (p < this->tot) subset.push_back(this->particles[p]);
        }
        return subset;
    }

    std::shared_ptr<Particle> random(){
        //https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
        //std::sample

        //return this->particles[(*distribution)(rand_gen)];
        return this->particles[Random::get_random(this->tot)];
    }


    template <typename T>
    void add(T com, T pos, double r, double rf, double q, double b, std::string name, bool image = false){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[this->tot]->com << com[0], com[1], com[2];
        this->particles[this->tot]->pos << pos[0], pos[1], pos[2];
        this->particles[this->tot]->qDisp = this->particles[this->tot]->pos - this->particles[this->tot]->com;
        this->particles[this->tot]->qDisp = this->particles[this->tot]->qDisp.stableNormalized() * b;
        //this->particles[this->tot]->pos = this->particles[this->tot]->com + this->particles[this->tot]->qDisp;

        //std::cout << this->particles[tot]->pos << " " << std::endl;
        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->name = name;


        if(!image){
            
            if(q > 0){
                this->cTot++;
            }
            else {
                this->aTot++;
            }

            this->pTot++;
        }
        else{
            this->iTot++;
        }

        this->tot++;
        assert(tot <= this->particles.size() && "tot is larger than particle vector size\n");
        assert(tot == this->pTot + this->iTot && "pTot + iTot != tot!\n");
    }

    template <typename T>
    void add(T com, double r, double rf, double q, double b, std::string name, bool image = false){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[this->tot]->com << com[0], com[1], com[2];
        this->particles[this->tot]->qDisp = Random::get_norm_vector();
        this->particles[this->tot]->qDisp = this->particles[this->tot]->qDisp.stableNormalized() * b;
        this->particles[this->tot]->pos = this->particles[this->tot]->com + this->particles[this->tot]->qDisp;
        //std::cout << this->particles[tot]->pos << " " << std::endl;
        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->name = name;


        if(!image){
            
            if(q > 0){
                this->cTot++;
            }
            else {
                this->aTot++;
            }

            this->pTot++;
        }
        else{
            this->iTot++;
        }

        this->tot++;
        assert(tot <= this->particles.size() && "tot is larger than particle vector size\n");
        assert(tot == this->pTot + this->iTot && "pTot + iTot != tot!\n");
    }


    void add(const std::shared_ptr<Particle>& p){

        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[tot]->com << p->com[0], p->com[1], p->com[2];
        this->particles[tot]->pos << p->pos[0], p->pos[1], p->pos[2];
        this->particles[tot]->qDisp << p->qDisp;

        this->particles[tot]->index = this->tot;
        this->particles[tot]->r = p->r;
        this->particles[tot]->rf = p->rf;
        this->particles[tot]->q = p->q;
        this->particles[tot]->b = p->b;
        this->particles[tot]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->pTot++;
        this->tot++;
    }


    void add(const std::shared_ptr<Particle>& p, int index){

        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }
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
        this->particles[index]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }      

        this->pTot++;
    }


    void remove(std::size_t index){
        //this->particles.erase(this->particles.begin() + i);
        //printf("Copying in remove\n");
        //std::copy(this->particles.begin() + i + 1, this->particles.begin() + this->tot, this->particles.begin() + i);
        //std::copy(this->particles.begin() + i, this->particles.begin() + this->tot, this->particles.begin() + i - 1);

        //this->particles.erase(this->particles.begin() + i);
        if(this->particles[index]->q > 0){
            this->cTot--;
        }
        else{
            this->aTot--; 
        }


        for(unsigned int i = index; i < this->tot - 1; i++){
            //printf("Remove: Moving particle %i to %i\n", i + 1, i);
            *(this->particles[i]) = *(this->particles[i + 1]);
            this->particles[i]->index = i;
        }

        this->pTot--;
        this->tot--;
    }



    void load(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<std::string> names){
        //assert correct sizes
        bool setNModel = false, setPModel = false;
        for(unsigned int i = 0; i < pos.size(); i++){

            this->add(com[i], pos[i], r[i], rf[i], charges[i], b[i], names[i]);
            //this->add(com[i], com[i], 2.5, rf[i], charges[i], 0.0, names[i]);
            
            if(!setPModel){
                if(charges[i] > 0){
                    pModel.q = charges[i];
                    pModel.b = b[i];
                    pModel.r = r[i];
                    pModel.rf = rf[i];
                    pModel.name = names[i];
                    setPModel = true;
                }
            }

            if(!setNModel){
                if(charges[i] < 0){
                    nModel.q = charges[i];
                    nModel.b = b[i];
                    nModel.r = r[i];
                    nModel.rf = rf[i];
                    nModel.name = names[i];
                    setNModel = true;
                }
            }

        }
        if(!setPModel || !setNModel){
            printf("Cation or anion model not set!\n");
            exit(1);
        }
        printf("Loaded %u particles, %u cations and %u anions.\n", this->tot, this->cTot, this->aTot);
    }



    void create(int pNum, int nNum, double p, double n, double rfp = 2.5, double rfn = 2.5, double rp = 2.5, double rn = 2.5, double bp = 0.0, double bn = 0.0){

        pModel.q = p;
        pModel.r = rp;
        pModel.rf = rfp;
        pModel.b = bp;
        pModel.name = "Na";

        nModel.q = n;
        nModel.r = rn;
        nModel.rf = rfn;
        nModel.b = bn;
        nModel.name = "Cl";

        Eigen::Vector3d com;
        for(int i = 0; i < pNum + nNum; i++){
            com = Random::get_vector();
            (i < pNum) ? this->add(com, rp, rfp, p, bp, "Na") : this->add(com, rn, rfn, n, bn, "Cl");
        }
        printf("\nCreated %i cations and %i anions\n", pNum, nNum);
    }
};