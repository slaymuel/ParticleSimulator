#pragma once

#include <Eigen/Dense>
#include <vector>
#include <memory>
#include "ran2_lib.h"
#include <iostream>
#include <fstream>
#include <random>

namespace Simulator{

namespace Random{    
    static constexpr double PI = 3.14159265359;

    inline namespace{
        std::random_device r;
        std::seed_seq ssq;
        std::default_random_engine rand_gen;
        std::unique_ptr< std::uniform_real_distribution<double> > real_dist = std::make_unique< std::uniform_real_distribution<double> >(0.0, 1.0);;
        std::vector< std::vector<double> > customVals;
    }

    void create_CDF(std::string filename);

    double get_random_from_distribution();
    double get_random();

    int get_random(int i);

    Eigen::Vector3d random_pos_box(double rf, std::vector<double> box);
    Eigen::Vector3d get_random_vector(double l);

    struct get_vector{
        operator Eigen::Vector3d();
        operator std::vector<double>();
    };
    struct get_norm_vector{
        operator Eigen::Vector3d();
        operator std::vector<double>();
    };
};

}