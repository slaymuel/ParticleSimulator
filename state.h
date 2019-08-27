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
        energy = (*energyFunc).all2all(this->particles);
        error = std::fabs((energy - cummulativeEnergy) / energy);

        if(this->particles.tot != _old->particles.tot){
            printf("_old state has %i particles and current state has %i.\n", _old->particles.tot, this->particles.tot);
            exit(1);
        }

        if(error > 1e-10){
            printf("\n\nEnergy drift is too large: %.12lf\n\n", error);
            exit(1);
        } 
    }

    void finalize(){
        _old = std::make_shared<State>();

        for(std::shared_ptr<Particle> p : this->particles.particles){
            _old->particles.add(p);
        }

        //Calculate the initial energy of the system
        this->energy = (*energyFunc).all2all(this->particles);
        this->cummulativeEnergy = this->energy;
    }

    void save(){
        for(auto i : this->movedParticles){
            if(i->index >= this->_old->particles.tot){
                this->_old->particles.add(i);
            }
            else{
                *(this->_old->particles.particles[i->index]) = *(this->particles.particles[i->index]);
            }
        }

        movedParticles.clear();
        cummulativeEnergy += this->dE;
    }


    void revert(){
        //Set moved partiles in current state equal to previous state
        //also need to set volume and maybe other properties
        for(auto i : this->movedParticles){
            if(this->particles.tot > _old->particles.tot){
                //this->particles.remove(this->particles.particles[this->particles.tot - 1]->index);
                this->particles.tot--;
                this->particles.pTot--;
            }
            else{
                *(this->particles.particles[i->index]) = *(_old->particles.particles[i->index]);
            }
        }
        movedParticles.clear();
    }


    //Get energy different between this and old state
    double get_energy_change(){
        for(auto p : movedParticles){
            if(!this->geo->is_inside(p->pos) || this->overlap(p->index)){
                //If moved outside box or overlap, return inf
                return std::numeric_limits<double>::infinity();
            }
        }

        double E1 = (*energyFunc)( this->_old->particles.get_subset(this->movedParticles), this->particles.particles );
        double E2 = (*energyFunc)( movedParticles, this->particles.particles );
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

    void add_particle(){
        //this->particles.add(this->geo->random_pos(), p->r, p->q, p->b, p->name);
    }

    void add_images(){
        for(int i = 0; i < this->particles.pTot; i++){
            this->particles.add(this->geo->mirror(this->particles.particles[i]->pos), this->particles.particles[i]->r, 
                                                 -this->particles.particles[i]->q,    this->particles.particles[i]->b, 
                                                  this->particles.particles[i]->name + "I", true);
        }
    }


    void equilibrate(){
        Eigen::Vector3d v;

        for(int i = 0; i < this->particles.pTot; i++){
            this->particles.particles[i]->pos = this->geo->random_pos();
        }
        
        int i = 0, overlaps = 1;
        Eigen::Vector3d oldPos;
        std::shared_ptr<Particle> p;

        //Move particles to prevent overlap
        while(overlaps > 0){
            p = this->particles.random();
            oldPos = p->pos;
            p->translate(10.0);

            if(this->overlap(p->index) || !this->geo->is_inside(p->pos)){
                p->pos = oldPos;
            }


            if(i % 50000 == 0){
                overlaps = this->get_overlaps();
                printf("Overlaps: %i, iteration: %i\r", overlaps, i);
                fflush(stdout);
            }
            i++;
            if(i > 1E9){
                i = 0;
            }
        }
        printf("\nEquilibration done\n\n");
    }



    bool overlap(std::size_t i){
        for(auto p : this->particles.particles){
            if(p->index == i) continue;

            if(this->geo->distance(p->pos, this->particles.particles[i]->pos) <= p->r + this->particles.particles[i]->r){
                return true;
            }
        }
        return false;
    }

    int get_overlaps(){
        int count = 0;
        for(auto p : this->particles.particles){
            (this->overlap(p->index)) ? count++ : 0;
        }
        return count;
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