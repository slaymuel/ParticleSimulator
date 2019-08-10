#pragma once

#include <Eigen/Dense>
#include <vector>
#include "particle.h"
#include <random>

class Particles{

    private:
    std::default_random_engine rand_gen;
    std::shared_ptr< std::uniform_int_distribution<int> > distribution;

    public:
    Eigen::MatrixXd positions;
    std::vector< std::shared_ptr<Particle> > particles;
    std::vector<int> movedParticles;
    std::vector<int> cations;
    std::vector<int> anions;
    
    template <typename T>
    void add(double x, double y, double z){
        this->positions.conservativeResize(this->positions.rows() + 1, 3);
        this->positions.row(this->positions.rows() - 1) << x, y, z;
        this->particles.push_back(std::make_shared<T>());
        this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles.back()->index = this->particles.size() - 1;
        //Update distribution for random generator
        distribution = std::make_shared< std::uniform_int_distribution<int> >(0, particles.size() - 1);
    }

    std::shared_ptr<Particle> get_random(){
        return particles[(*distribution)(rand_gen)];
    }

    Eigen::MatrixXd get_subset(int sr, int fr){
        return this->positions.block(sr, 0, fr, 3);
    }
};