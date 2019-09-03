#pragma once

#include "particle.h"
#include "particles.h"

class EnergyBase{
    public:
    virtual ~EnergyBase(){};
    virtual double all2all(Particles& particles) = 0;
    virtual double i2all(std::shared_ptr<Particle> p, Particles& particles) = 0;
    virtual double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) = 0;
    virtual double operator()(std::vector< int >&& p, Particles& particles) = 0;
    virtual double operator()(std::vector< int >& p, Particles& particles) = 0;
};


template <typename E>
class Energy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    double all2all(Particles& particles){
        double e = 0.0;
        for(int i = 0; i < particles.tot; i++){
            for(int j = i + 1; j < particles.tot; j++){
                e += i2i(particles.particles[i], particles.particles[j]);
            } 
        }
        return e * constants::lB;
    }

    double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = 0.0;
        for(int i = 0; i < particles.tot; i++){
            if(p->index == particles.particles[i]->index) continue;
            e += i2i(p, particles.particles[i]);
        }
        return e * constants::lB;
    }

    double operator()(std::vector< int >&& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e;
    }

    double operator()(std::vector< int >& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e;
    }

    inline double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        return energy_func(p1, p2);
    }
};