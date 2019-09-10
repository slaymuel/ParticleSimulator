#pragma once

#include "particle.h"
#include "particles.h"
#include "geometry.h"


class EnergyBase{
    protected:
    Geometry *geo;

    public:
    void set_geo(Geometry* geo){
        this->geo = geo;
    }

    virtual ~EnergyBase(){};
    virtual double all2all(Particles& particles) = 0;
    virtual double i2all(std::shared_ptr<Particle> p, Particles& particles) = 0;
    virtual double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) = 0;
    virtual double operator()(std::vector< int >&& p, Particles& particles) = 0;
    virtual double operator()(std::vector< int >& p, Particles& particles) = 0;
    virtual void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new) = 0;
    virtual void initialize(Particles& particles) = 0;
};


template <typename E>
class PairEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    double all2all(Particles& particles){
        double e = 0.0;


        #pragma omp parallel for reduction(+:e) schedule(guided, 50) if(particles.tot >= 500)
        for(int i = 0; i < particles.tot; i++){
            for(int j = i + 1; j < particles.tot; j++){
                e += i2i(particles.particles[i], particles.particles[j]);
            } 
        }

        return e * constants::lB;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = 0.0;

        #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
        for (int i = 0; i < particles.tot; i++){
            if (p->index == particles.particles[i]->index) continue;
            e += i2i(p, particles.particles[i]);
        }

        return e;
    }

    double operator()(std::vector< int >&& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }

        return e * constants::lB;
    }

    double operator()(std::vector< int >& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e * constants::lB;
    }

    inline double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        return energy_func(p1, p2, geo->distance(p1->pos, p2->pos));
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    virtual void initialize(Particles& particles){}
};




template <typename E>
class ExtEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:
    ExtEnergy(double x, double y, double z){
        energy_func.set_box(x, y, z);
    }

    double i2all(std::shared_ptr<Particle> p, Particles& particles){ return 0.0; }
    double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){ return 0.0; }

    double all2all(Particles& particles){
        return energy_func() * constants::lB;
    }

    double operator()(std::vector< int >&& p, Particles& particles){
        return energy_func() * constants::lB;
    }

    double operator()(std::vector< int >& p, Particles& particles){
        return energy_func() * constants::lB;
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){
        energy_func.update(_old, _new);
    }
    void initialize(Particles& particles){
        energy_func.initialize(particles);
    }
};