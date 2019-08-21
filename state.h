#pragma once

#include <vector>
#include "particles.h"
#include <limits>
//#include <math.h>
#include "geometry.h"
#include "energy.h"
#include "potentials.h"

class State{
    private:

    std::shared_ptr<State> _old;

    public:

    double energy, cummulativeEnergy, dE, error;
    Particles particles;
    std::vector< std::shared_ptr<Particle> > movedParticles;    //Particles that has moved from previous state
    Geometry *geo;
    std::shared_ptr<EnergyBase> energyFunc;
    

    void control(){
        energy = (*energyFunc).all2all(this->particles.particles);
        error = std::fabs((energy - cummulativeEnergy) / energy);

        if(error > 1e-10){
            printf("\n\nEnergy drift is too large!\n\n");
            exit(1);
        } 
    }

    void finalize(){
        _old = std::make_shared<State>();

        for(std::shared_ptr<Particle> p : this->particles.particles){
            _old->particles.particles.push_back(std::make_shared<Particle>());
            *(_old->particles.particles.back()) = *p;
        }

        //Calculate the initial energy of the system
        this->energy = (*energyFunc).all2all(this->particles.particles);
        this->cummulativeEnergy = this->energy;
    }

    void save(){
        for(auto i : this->movedParticles){
            *(_old->particles.particles[i->index]) = *(this->particles.particles[i->index]);
        }
        movedParticles.clear();
        cummulativeEnergy += this->dE;
    }


    void revert(){
        //Set moved partiles in current state equal to previous state
        //also need to set volume and maybe other properties
        for(auto i : this->movedParticles){
            *(this->particles.particles[i->index]) = *(_old->particles.particles[i->index]);
        }
        movedParticles.clear();
    }


    //Get energy different between this and old state
    double get_energy_change(){
        for(auto p : movedParticles){
            if(!geo->is_inside(this->particles.particles[p->index]->pos) || this->overlap(p->index)){
                //If moved outside box or overlap, return inf
                return std::numeric_limits<double>::infinity();
            }
        }

        double E1 = (*energyFunc)( this->_old->particles.get_subset(this->movedParticles), this->particles.particles );
        double E2 = (*energyFunc)(movedParticles, this->particles.particles);
        this->dE = E2 - E1;

        return dE;
    }


    //Called when a move is accepted - set movedParticles
    void move_callback(std::vector< std::shared_ptr<Particle> > ps){   
        //this->movedParticles.insert(std::end(movedParticles), std::begin(ps), std::end(ps));
        std::for_each(std::begin(ps), std::end(ps), [this](std::shared_ptr<Particle> p){ 
                                                            this->movedParticles.push_back(p); 
                                                            });
    }


    void load_state(std::vector<double> pos){
        for(auto p : pos){
            std::cout << p << std::endl;
        }
    }


    void equilibrate(){
        std::vector<double> v;
        for(auto p : this->particles.particles){
            v = Random::get_vector();
            p->pos[0] = this->geo->dh[0] * (v[0] * 2.0 - 1);
            p->pos[1] = this->geo->dh[1] * (v[1] * 2.0 - 1);
            p->pos[2] = this->geo->dh[2] * (v[2] * 2.0 - 1);
        }
    }


    bool overlap(std::size_t i){
        for(auto p : this->particles.particles){
            if(p->index == i) continue;

            if(geo->distance(p->pos, this->particles.particles[i]->pos) <= p->r + this->particles.particles[i]->r){
                return true;
            }
        }
        return false;
    }


    void set_geometry(int type){

        switch (type){
            default:
                printf("Creating Cuboid box\n");
                this->geo = new Cuboid<true, true, true>(50.0, 50.0, 50.0);
                break;
            case 1:
                this->geo = new Sphere();
                break;
        }
    }
    


    void set_energy(int type){
        switch (type){
            default:
                printf("Adding Coulomb potential\n");
                this->energyFunc = std::make_shared< Energy<Coulomb> >();
                break;
        }
        
    }
};