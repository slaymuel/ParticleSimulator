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
        printf("Saving state\n");
        for(auto i : this->movedParticles){
            *(_old->particles.particles[i->index]) = *(this->particles.particles[i->index]);
        }
        movedParticles.clear();
        cummulativeEnergy += this->dE;
    }


    void revert(){
        printf("Reverting state\n");
        //Set moved partiles in current state equal to previous state
        //also need to set volume and maybe other properties
        for(auto i : this->movedParticles){
            printf("new %lf %lf %lf, old %lf %lf %lf\n", this->particles.particles[i->index]->pos[0], this->particles.particles[i->index]->pos[1], this->particles.particles[i->index]->pos[2],
                                                    _old->particles.particles[i->index]->pos[0], _old->particles.particles[i->index]->pos[1], _old->particles.particles[i->index]->pos[2]);
            *(this->particles.particles[i->index]) = *(_old->particles.particles[i->index]);
        }
        movedParticles.clear();
    }

    //Get energy different between this and old state
    double get_energy_change(){ 
        double E1 = 0.0, E2 = 0.0;

        for(auto p : movedParticles){
            if(!geo->is_inside(this->particles.particles[p->index]->pos)){
                //If moved outside box or overlap, return inf
                return std::numeric_limits<double>::infinity();
            }
        }
        E1 += (*energyFunc)(_old->particles.get_subset(this->movedParticles), this->particles.particles);
        E2 += (*energyFunc)(movedParticles, this->particles.particles);
        this->dE = E2 - E1;

        printf("Energies = new %lf, old %lf\n", E2, E1);

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