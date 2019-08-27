#pragma once

#include "particle.h"
#include "particles.h"

class EnergyBase{
    public:
    virtual ~EnergyBase(){};
    virtual double all2all(Particles& particles) = 0;
    virtual double i2all(std::shared_ptr<Particle> p, std::vector< std::shared_ptr<Particle> >& particles) = 0;
    virtual double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) = 0;
    virtual double operator()(std::vector< std::shared_ptr<Particle> >&& p, std::vector< std::shared_ptr<Particle> >& particles) = 0;
    virtual double operator()(std::vector< std::shared_ptr<Particle> >& p, std::vector< std::shared_ptr<Particle> >& particles) = 0;
};


template <typename E>
class Energy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    double all2all(Particles& particles){
        double energy = 0;
        for(int i = 0; i < particles.tot; i++){
            for(int j = i + 1; j < particles.tot; j++){
                energy += i2i(particles.particles[i], particles.particles[j]);
            } 
        }
        return energy * constants::lB;
    }

    double i2all(std::shared_ptr<Particle> p, std::vector< std::shared_ptr<Particle> >& particles){
        double energy = 0;
        for(auto particle : particles){
            if(p->index == particle->index) continue;
            energy += i2i(p, particle);
        }
        return energy * constants::lB;
    }

    double operator()(std::vector< std::shared_ptr<Particle> >&& p, std::vector< std::shared_ptr<Particle> >& particles){
        double e = 0;
        for(auto s : p){
            e += i2all(s, particles);
        }
        return e;
    }

    double operator()(std::vector< std::shared_ptr<Particle> >& p, std::vector< std::shared_ptr<Particle> >& particles){
        return i2all(p[0], particles);
    }

    inline double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        return energy_func(p1, p2);
    }
};