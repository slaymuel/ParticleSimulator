#pragma once
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include "constants.h"
#include "ran2_lib.h"
#include <iostream>
#include <fstream>

class Random{
    private:
    static std::random_device r;
    static std::seed_seq ssq;
    static std::default_random_engine rand_gen;
    static std::unique_ptr< std::uniform_real_distribution<double> > real_dist;
    static std::vector< std::vector<double> > customVals;

    public:

    static void initialize(){

        real_dist = std::make_unique< std::uniform_real_distribution<double> >(0.0, 1.0);
        //create_CDF("pqDist_ideal.txt");
    }

    static void create_CDF(std::string filename){
        
        
        std::cout << "Creating CDF" << std::endl;
        std::ifstream f(filename);
        std::string line;
        double integral = 0;
        
        while (std::getline(f, line)) {
            double x, y;
            
            std::istringstream ss(line);
            ss >> x >> y;
            std::cout << x << " " << y << std::endl;

            if(y > 1e-10){
                customVals.push_back( std::vector<double>({x, y}) );
                integral += y;
            }
        }
        double dx = customVals[1][0] - customVals[0][0];
        std::cout << "dx is: " << dx << std::endl;
        //integral *= dx;
        for(std::vector<double>& val : customVals){
            val[1] /= integral;
        }

        for(int i = 1; i < customVals.size(); i++){
            customVals[i][1] += customVals[i - 1][1];
            printf("vals %lf %lf\n", customVals[i][0], customVals[i][1]);
        }
    }

    inline static double get_random_from_distribution(){
        double r = get_random();

        for(int i = 0; i < customVals.size(); i++){
            if(r < customVals[i][1]){
                return customVals[i][0];
            }
        }
        return -1.0;
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

    inline static Eigen::Vector3d get_random_vector(double l){
        Eigen::Vector3d v(get_random()*2.0*l - l, get_random()*2.0*l - l, get_random()*2.0*l - l);
        if(l < 1e-10){
            return v;
        }
        while(v.norm() > l){
            v << get_random()*2.0*l - l, get_random()*2.0*l - l, get_random()*2.0*l - l;
        }
        return v;
    }

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
        
        inline operator std::vector<double>(){
            
            double phi = get_random() * 2.0 * constants::PI;
            double z = get_random() * 2.0 - 1.0;
            std::vector<double> v = {std::sqrt(1.0 - z*z) * std::cos(phi), std::sqrt(1.0 - z*z) * std::sin(phi), z};
            
            return v;
        }
    };
};


std::unique_ptr< std::uniform_real_distribution<double> > Random::real_dist;
std::random_device Random::r;
std::seed_seq Random::ssq{r()};
std::default_random_engine Random::rand_gen(ssq);
std::vector< std::vector<double> > Random::customVals;