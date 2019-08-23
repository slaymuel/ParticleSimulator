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
    std::vector< std::shared_ptr<Particle> > particles;
    std::vector<int> movedParticles;
    std::vector<int> cations;
    std::vector<int> anions;
    int pTot;

    //Eigen::MatrixXd get_subset(int sr, int fr){
    //    return this->positions.block(sr, 0, fr, 3);
    //}

    Particles(){
        this->pTot = 0;
    }

    std::vector< std::shared_ptr<Particle> > get_subset(std::vector< std::shared_ptr<Particle> >& ps){
        std::vector< std::shared_ptr<Particle> > subset;
        for(std::shared_ptr<Particle> p : ps){
            subset.push_back(particles[p->index]);
        }
        return subset;
    }

    std::shared_ptr<Particle> get_random(){
        //https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
        //std::sample

        return this->particles[(*distribution)(rand_gen)];
    }



    void add(std::vector<double> pos, double r, double q, double b, std::string name){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];
        //Create particle
        if(pTot >= this->particles.size()){
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[pTot]->pos  << pos[0], pos[1], pos[2];
        this->particles[pTot]->index = this->particles.size() - 1;
        this->particles[pTot]->r = r;
        this->particles[pTot]->q = q;
        this->particles[pTot]->b = b;
        this->particles[pTot]->name = name;
        //Update distribution for random generator
        this->distribution = std::make_shared< std::uniform_int_distribution<int> >(0, particles.size() - 1);
        this->pTot++;
    }


    void add(std::shared_ptr<Particle> p){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];
        //Create particle
        if(pTot >= this->particles.size()){
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[pTot]->pos  << p->pos[0], p->pos[1], p->pos[2];
        this->particles[pTot]->index = this->particles.size() - 1;
        this->particles[pTot]->r = p->r;
        this->particles[pTot]->q = p->q;
        this->particles[pTot]->b = p->b;
        this->particles[pTot]->name = p->name;
        //Update distribution for random generator
        this->distribution = std::make_shared< std::uniform_int_distribution<int> >(0, particles.size() - 1);
        this->pTot++;
    }


    void remove(std::size_t i){
        //this->particles.erase(this->particles.begin() + i);
        std::copy(this->particles.begin() + i, this->particles.begin() + this->pTot, this->particles.begin() + i - 1);
        for(;i < this->particles.size(); i++){
            this->particles[i]->index--;
        }
        this->pTot--;
    }



    void load(std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> b, std::vector<std::string> names){
        //assert correct sizes
        for(int i = 0; i < pos.size(); i++){
            this->add(pos[i], 2.5, charges[i], b[i], names[i]);
        }
    }



    void create(int pNum, int nNum){
        int tot = pNum + nNum;
        std::vector<double> v;

        for(int i = 0; i < tot; i++){
            v = Random::get_vector();
            (i % 2 == 0) ? this->add(v, 2.5, 1.0, 0.0, "Na") : this->add(v, 2.5, -1.0, 0.0, "Cl") ;
        }
    }


    //Write xyz file
    void to_xyz(std::string fileName){
        std::ofstream f (fileName);
        if (f.is_open())
        {
            f << this->particles.size() << "\n\n";
            for(auto p : this->particles){

                f << std::fixed << std::setprecision(3) << p->name << " " <<  p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << "\n";
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