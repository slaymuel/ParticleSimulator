#pragma once

#include <vector>
#include "particles.h"
#include <limits>
//#include <math.h>
#include "geometry.h"
#include "energy.h"
#include "potentials.h"

class State{
    public:

    double energy;
    Particles particles;
    std::vector< std::shared_ptr<Particle> > movedParticles;    //Particles that has moved from previous state
    Geometry *geo;
    std::shared_ptr<EnergyBase> energyFunc;
    
    void save(){
        movedParticles.clear();
    }


    void revert(State& old){
        //Set moved partiles in current state equal to previous state
        //also need to set volume and maybe other properties
        for(auto i : this->movedParticles){
            this->particles.particles[i->index] = old.particles.particles[i->index];
        }
        movedParticles.clear();
    }


    double get_energy_change(State old){ //get energy different between this and old state
        double dE = 0.0;

        for(auto p : movedParticles){
            if(!geo->is_inside(this->particles.particles[p->index]->pos)){
                return std::numeric_limits<double>::infinity();
            }
            dE += (*energyFunc)(movedParticles, this->particles.particles);
        }
        //If moved outside box or overlap, return inf
        //dE = std::numeric_limits<double>::infinity();
        printf("Energy = %lf\n", dE);
        return dE;
    }


    void move_callback(std::vector< std::shared_ptr<Particle> > ps){   //Called when a move is accepted
        //this->movedParticles.insert(std::end(movedParticles), std::begin(indices), std::end(indices));
        std::for_each(std::begin(ps), std::end(ps), [this](std::shared_ptr<Particle> p){ 
                                                            this->movedParticles.push_back(p); 
                                                            });
        //for(auto p : ps){
        //    this->movedParticles.push_back(p->index);
            //particles.particles[p->index] = p; //dont need to do this since state.particle is already moved
        //}
    }


    void load_state(std::vector<double> pos){
        for(auto p : pos){
            std::cout << p << std::endl;
        }
    }


    void load_particles(std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> b){
        //assert correct sizes
        for(int i = 0; i < pos.size(); i++){
            particles.add<Particle>(pos[i], charges[i], b[i]);
        }
    }

    void set_geometry(int type){

        switch (type){
            default:
                printf("Creating Cuboid box\n");
                this->geo = new Cuboid(10.0, 10.0, 10.0);
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