#pragma once

#include <Eigen/Dense>
#include <vector>
#include "particle.h"
#include <chrono>

class Particles{

    private:
    std::default_random_engine rand_gen;
    std::unique_ptr< std::uniform_int_distribution<int> > distribution;

    public:
    Eigen::MatrixXd positions;
    std::vector< std::shared_ptr<Particle> > particles;
    std::vector<int> movedParticles;
    std::vector<int> cations;
    std::vector<int> anions;
    
    template <typename T>
    void add(std::vector<double> pos, double q, double b){
        //Resize positions
        this->positions.conservativeResize(this->positions.rows() + 1, 3);
        this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];
        //Create particle
        this->particles.push_back(std::make_shared<T>());
        this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles.back()->index = this->particles.size() - 1;
        this->particles.back()->q = q;
        this->particles.back()->b = b;

        //Update distribution for random generator
        distribution = std::make_unique< std::uniform_int_distribution<int> >(0, particles.size() - 1);
    }

    std::shared_ptr<Particle> get_random(){
        //https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
        //std::sample

        return particles[(*distribution)(rand_gen)];
    }

    //Eigen::MatrixXd get_subset(int sr, int fr){
    //    return this->positions.block(sr, 0, fr, 3);
    //}
    std::vector< std::shared_ptr<Particle> > get_subset(std::vector< std::shared_ptr<Particle> >& ps){
        std::vector< std::shared_ptr<Particle> > subset;
        for(std::shared_ptr<Particle> p : ps){
            subset.push_back(particles[p->index]);
        }
        return subset;
    }

    bool overlap(std::size_t i){
        for(auto p : this->particles){
            if(p->index == i) continue;

            if(p->distance(this->particles[i]) <= p->r + this->particles[i]->r){
                return true;
            }
        }
        return false;
    }
};