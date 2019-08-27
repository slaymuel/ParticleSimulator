#pragma once

#include <Eigen/Dense>
#include <vector>
#include "particle.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

class Particles{

    private:
    std::random_device r;
    std::seed_seq ssq{r()};
    std::default_random_engine rand_gen{ssq};
    std::shared_ptr< std::uniform_int_distribution<int> > distribution;

    public:
    //Eigen::MatrixXd positions;
    std::vector< std::shared_ptr<Particle> > particles, cations, anions;
    std::vector<int> movedParticles;
    int pTot = 0, cTot = 0, aTot = 0, iTot = 0, tot = 0;

    //Eigen::MatrixXd get_subset(int sr, int fr){
    //    return this->positions.block(sr, 0, fr, 3);
    //}

    Particles(){}

    std::vector< std::shared_ptr<Particle> > get_subset(std::vector< std::shared_ptr<Particle> >& ps){
        std::vector< std::shared_ptr<Particle> > subset;
        for(std::shared_ptr<Particle> p : ps){
            if(p->index < this->tot) subset.push_back(this->particles[p->index]);
        }
        return subset;
    }

    std::shared_ptr<Particle> random(){
        //https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
        //std::sample

        return this->particles[(*distribution)(rand_gen)];
    }


    template <typename T>
    void add(T pos, double r, double q, double b, std::string name, bool image = false){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        if(tot == this->particles.size()){
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[tot]->pos  << pos[0], pos[1], pos[2];
        this->particles[tot]->index = this->tot;
        this->particles[tot]->r = r;
        this->particles[tot]->q = q;
        this->particles[tot]->b = b;
        this->particles[tot]->name = name;
        //Update distribution for random generator

        if(!image){
            std::shared_ptr<Particle> temp = this->particles[pTot];
            this->pTot++;
            if(q > 0){
                this->cations.push_back(temp);
                this->cTot++;
            }
            else {
                this->anions.push_back(temp);
                this->aTot++;
            }

            this->distribution = std::make_shared< std::uniform_int_distribution<int> >(0, this->pTot - 1);
        }
        else{
            this->iTot++;
        }

        tot++;
        assert(tot <= this->particles.size() && "tot is larger than particle vector size\n");
        assert(tot == pTot + iTot && "pTot + iTot != tot!\n");
    }


    void add(std::shared_ptr<Particle> p){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];
        //Create particle
        if(tot == this->particles.size()){
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[tot]->pos  << p->pos[0], p->pos[1], p->pos[2];
        this->particles[tot]->index = this->tot;
        this->particles[tot]->r = p->r;
        this->particles[tot]->q = p->q;
        this->particles[tot]->b = p->b;
        this->particles[tot]->name = p->name;
        this->pTot++;
        //Update distribution for random generator
        this->distribution = std::make_shared< std::uniform_int_distribution<int> >(0, this->pTot - 1);

        this->tot++;
    }


    void remove(std::size_t i){
        //this->particles.erase(this->particles.begin() + i);
        std::copy(this->particles.begin() + i, this->particles.begin() + this->tot, this->particles.begin() + i - 1);
        for(;i < this->particles.size(); i++){
            this->particles[i]->index--;
        }
        this->pTot--;
        this->tot--;
    }



    void load(std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> b, std::vector<std::string> names){
        //assert correct sizes
        for(int i = 0; i < pos.size(); i++){
            this->add(pos[i], 2.5, charges[i], b[i], names[i]);
        }
    }



    void create(int pNum, int nNum){
        std::vector<double> v;

        for(int i = 0; i < pNum + nNum; i++){
            v = Random::get_vector();
            (i < pNum == 0) ? this->add(v, 2.5, 1.0, 0.0, "Na") : this->add(v, 2.5, -1.0, 0.0, "Cl") ;
        }
    }


    //Write xyz file
    void to_xyz(std::string fileName){
        std::ofstream f (fileName);
        if (f.is_open())
        {
            f << this->tot << "\n\n";
            for(int i = 0; i < this->tot; i++){

                f << std::fixed << std::setprecision(3) << this->particles[i]->name << " " <<  this->particles[i]->pos[0] << " " << this->particles[i]->pos[1] << " " << this->particles[i]->pos[2] << "\n";
            }
            f << "10 10 10" << "\n";
            f.close();
        }
        else std::cout << "Unable to open file";
    }

    //write checkpoint file
    void to_cpt(){

    }
};