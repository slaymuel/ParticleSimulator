#pragma once

#include <vector>
#include "particles.h"
#include <limits>
//#include <math.h>

class State{
    public:

    
    double energy;
    Particles particles;
    std::vector<int> movedParticles;    //Particles that has moved from previous state
    //Geometry

    void save(){
        for(auto i : this->movedParticles){
            //Update current state
            //Dont need to do anythin since 
            //current state is already updated
        }
        movedParticles.clear();
    }

    void revert(State& old){
        //Set moved partiles in current state equal to previous state
        //also need to set volume and maybe other properties
        for(auto i : this->movedParticles){
            this->particles.particles[i] = old.particles.particles[i];
        }
        movedParticles.clear();
    }

    double get_energy_change(State old){ //get energy different between this and old state
        double dE = 0.0;

        for(int particle : movedParticles){
            //dE += get_energy(particle);
        }
        //If moved outside box or overlap, return inf
        //dE = std::numeric_limits<double>::infinity();
        return dE;
    }

    void move_callback(std::vector< std::shared_ptr<Particle> > ps){   //Called when a move is accepted
        //this->movedParticles.insert(std::end(movedParticles), std::begin(indices), std::end(indices));
        std::for_each(std::begin(ps), std::end(ps), [this](std::shared_ptr<Particle> p){ 
                                                            this->movedParticles.push_back(p->index); 
                                                            });
        //for(auto p : ps){
        //    this->movedParticles.push_back(p->index);
            //particles.particles[p->index] = p; //dont need to do this since state.particle is already moved
        //}
    }
};