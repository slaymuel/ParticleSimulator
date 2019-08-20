#pragma once
#include <Eigen/Dense>

class Random{
    private:

    static std::default_random_engine rand_gen;
    static std::shared_ptr< std::uniform_real_distribution<double> > real_dist;
    
    public:

    static void initialize(){
        real_dist = std::make_shared< std::uniform_real_distribution<double> >(0.0, 1.0);
    }

    static double get_random(){
        return (*real_dist)(rand_gen);
    }

    static Eigen::Vector3d get_random_vector(){
        Eigen::Vector3d a((double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen));
        return a.normalized();
    }
};

std::default_random_engine Random::rand_gen;
std::shared_ptr< std::uniform_real_distribution<double> > Random::real_dist;