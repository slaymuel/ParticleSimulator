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


        #pragma omp parallel for reduction(+:e) schedule(guided, 100) if(particles.tot >= 500)
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
    void initialize(Particles& particles){}
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




template <typename E>
class ImgEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double CC = 0.0, CpC = 0.0, self = 0.0;
        Eigen::Vector3d temp;

        // CC
        for (int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;

            CC += energy_func(p->q, particles[i]->q, this->geo->distance(p->pos, particles[i]->pos));
        }


        //  C'C
        temp = p->pos;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        for (int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;

            CpC += energy_func(-p->q, particles[i]->q, this->geo->distance(temp, particles[i]->pos));
        }
        // => CC == C'C' and C'C == CC'

        // Self term  
        self = energy_func(p->q, -p->q, this->geo->distance(p->pos, temp));
        return CC + CpC + 0.5 * self;
    }


    inline double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        return energy_func(p1->q, p2->q, this->geo->distance(p1->pos, p2->pos));
    }


    double all2all(Particles& particles){
        double CC = 0.0, CpC = 0.0;
        Eigen::Vector3d temp;
        Eigen::Vector3d temp2;

        // CC
        //#pragma omp parallel for schedule(guided, 100) if(particles.tot >= 500)
        for(int i = 0; i < particles.tot; i++){
            for(int j = i + 1; j < particles.tot; j++){
                CC += energy_func(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
            } 
        }

        //C'C
        //#pragma omp parallel for schedule(dynamic, 100) reduction(+:CpC) private(temp)
        for(int i = 0; i < particles.tot; i++){
            temp = particles[i]->pos;
            temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 

            for(int j = 0; j < particles.tot; j++){
                CpC += energy_func(-particles[i]->q, particles[j]->q, this->geo->distance(temp, particles[j]->pos));
            } 
        }

        return (CC + 0.5 * CpC) * constants::lB;
    }


    double operator()(std::vector< int >&& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles[s], particles);
        }
        return e * constants::lB;
    }


    double operator()(std::vector< int >& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles[s], particles);
        }
        return e * constants::lB;
    }


    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){
        //energy_func.update(_old, _new);
    }


    void initialize(Particles& particles){
        //energy_func.initialize(particles);
    }


};