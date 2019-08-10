
#include <iostream>
#include "particles.h"
#include "state.h"
#include "particle.h"
#include <vector>
#include "move.h"
#include "simulator.h"

int main(){
    /*std::vector<Particle*> ps;
    Particles particles;
    for(int i = 0; i < 10; i++){
        particles.add(i,i,i);
    }*/
    //Move* move = new Translate();
    //ps.push_back(new RPM());
    //ps.push_back(new ARPM());
    //ps[0]->translate();
    //ps[1]->translate();
    //trans.operator()<decltype(ps[1])>(ps[0]);

    State currentState;
    Simulator* sim = new Simulator(10, 10);
    sim->run();
    //std::function<void(std::vector<int>)> move_callback = std::bind(&State::move_callback, &currentState, std::placeholders::_1);
    //std::function<void(std::vector<int>)> move_callback = [currentState](std::vector<int> indices) { currentState.move_callback(indices); }
    /*std::cout << "size before: " << currentState.movedParticles.size() << std::endl;
    (*move)(ps[1], move_callback);
    std::cout << "size after: " <<  currentState.movedParticles.size() << std::endl;


    std::cout << particles.positions << std::endl;
    std::cout << particles.get_subset(1, 3).transpose() << std::endl;*/


    //Simulator sim;
    //sim.setup(params);
    //sim.run();
}