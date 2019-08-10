#pragma once

#include "state.h"
#include "move.h"
#include <functional>

class Simulator{
    private:
    int macroSteps;
    int microSteps;
    std::vector<Particle*> ps;
    State previousState;
    State currentState;
    //Particles *particles = new Particles();
    
    public:
    Simulator(int macro, int micro) : macroSteps(macro), microSteps(micro){}
    //energy
    //moves
    
    std::vector<Move*> moves;
    
    /* State callback after move */
    std::function< void(std::vector< std::shared_ptr<Particle> >) > move_callback 
                = std::bind(&State::move_callback, &currentState, std::placeholders::_1);
    
    void run(){
        currentState.particles.add<RPM>(1,1,1);
        currentState.particles.add<ARPM>(1,1,1);
        ps.push_back(new RPM());
        ps.push_back(new ARPM());
        moves.push_back(new Translate());
        moves.push_back(new Rotate());

        for(int macro = 0; macro < macroSteps; macro++){
            for(int micro = 0; micro < microSteps; micro++){
                for(auto move : moves){
                    std::cout << "moved before: " << currentState.movedParticles.size() << std::endl;

                    //Move  
                    //Move should check if particle is part of molecule
                    (*move)(currentState.particles.get_random(), move_callback); // Two virtual calls
                    //state->update();
                    std::cout << "moved after: " <<  currentState.movedParticles.back() << std::endl;
                                                            

                    if(move->accept( currentState.get_energy_change(previousState) )){
                        currentState.save();
                    }
                    else{
                        currentState.revert(previousState);
                    }
                    //should also be able to
                    //state.get_energy(subset_of_particles);
                    //energy.get_energy(subset of particles)
                }
            }
            //1. Lista/vektor med olika input som de olika samplingsmetoderna behöver
            //2. sampler kan på något sätt efterfråga input, text genom att sätta en variabel
            //   Sen kan simulator ha en map och leta på den variabeln
            //sampler.sample(state);
            
        }
    }
};