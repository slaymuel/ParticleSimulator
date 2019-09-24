#pragma once
#include <Eigen/Dense>
#include <vector>
#include <memory>

class Random{
    private:
    static std::random_device r;
    static std::seed_seq ssq;
    static std::default_random_engine rand_gen;
    static std::unique_ptr< std::uniform_real_distribution<double> > real_dist;
    
    public:

    static void initialize(){

        real_dist = std::make_unique< std::uniform_real_distribution<double> >(0.0, 1.0);
    }

    inline static double get_random(){
        return (*real_dist)(rand_gen);
    }

    inline static int get_random(int i){
        return i * (*real_dist)(rand_gen);
    }

    inline static Eigen::Vector3d get_random_vector(){
        Eigen::Vector3d a((double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen));
        return a.normalized();
    }

    struct get_vector{
        inline operator Eigen::Vector3d(){
            Eigen::Vector3d v((double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen));
            return v;
        }
        inline operator std::vector<double>(){
            std::vector<double> v = { (double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen) };
            return v;
        }
    };
};


std::unique_ptr< std::uniform_real_distribution<double> > Random::real_dist;
std::random_device Random::r;
std::seed_seq Random::ssq{r()};
std::default_random_engine Random::rand_gen(ssq);