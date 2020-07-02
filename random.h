#pragma once
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include "constants.h"
#include "ran2_lib.h"

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
        //return (*real_dist)(rand_gen);
        return ran2::get_random();
    }

    inline static int get_random(int i){
        return i * get_random();
    }

    inline static Eigen::Vector3d random_pos_box(double rf, std::vector<double> box){
        Eigen::Vector3d v;
        v = get_vector();
        v << box[0] * v[0], box[1] * v[1], (box[2] - rf) * v[2];
        return v;
    }

    //inline static Eigen::Vector3d get_random_vector(){
    //    Eigen::Vector3d a((double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen), (double)(*real_dist)(rand_gen));
    //    return a;
    //}

    struct get_vector{
        inline operator Eigen::Vector3d(){
            Eigen::Vector3d v(get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0);
            return v;
        }
        inline operator std::vector<double>(){
            std::vector<double> v = {get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0};
            return v;
        }
    };

    struct get_norm_vector{
        inline operator Eigen::Vector3d(){
            double phi = get_random() * 2.0 * constants::PI;
            double z = get_random() * 2.0 - 1.0;
            Eigen::Vector3d v(std::sqrt(1.0 - z*z) * std::cos(phi), std::sqrt(1.0 - z*z) * std::sin(phi), z);
            return v;
        }
    };
};


std::unique_ptr< std::uniform_real_distribution<double> > Random::real_dist;
std::random_device Random::r;
std::seed_seq Random::ssq{r()};
std::default_random_engine Random::rand_gen(ssq);