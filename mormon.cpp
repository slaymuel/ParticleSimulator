#include <random>
#include "random.h"
#include <assert.h>
#include "constants.h"
#include <iostream>
#include "state.h"
#include "particle.h"
#include "move.h"
#include <functional>

#ifdef PY11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
#endif



class Simulator{
    private:
    std::vector<Particle*> ps;
    
    public:
    State state;
    Simulator(double Dielec, double T){
        //Set some constants
        constants::D = Dielec;
        constants::T = T; 
        constants::lB = constants::C * (1.0 / (Dielec * T));
        Random::initialize();
    }

    //Energy<Coulomb> energy;
    //moves
    
    std::vector<Move*> moves;
    
    /* State callback after move */
    std::function< void(std::vector< std::shared_ptr<Particle> >) > move_callback 
                = std::bind(&State::move_callback, &state, std::placeholders::_1);


    void run(int macroSteps, int microSteps){
        std::cout << "Running simulation at: " << constants::T << "K with: " << state.particles.particles.size() 
                                                                << " particles" << std::endl;

        moves.push_back(new Translate(0.2));
        //moves.push_back(new Rotate());

        for(int macro = 0; macro < macroSteps; macro++){
            for(int micro = 0; micro < microSteps; micro++){
                for(auto move : moves){

                    //Move should check if particle is part of molecule
                    (*move)(state.particles.random(), move_callback); // Two virtual calls
                    
                    if(move->accept( state.get_energy_change() )){
                        state.save();
                    }
                    else{
                        state.revert();
                    }

                    //should also be able to
                    //state.get_energy(subset_of_particles);
                }
            }

            /*                                HALF TIME                                  */

            //Check energy drift etc
            state.control();

            //Print progress
            std::cout << "\nIteration: " << macro * microSteps + microSteps << std::endl;
            printf("Acceptance: %.1lf%%\n", (double)moves[0]->accepted / Move::totalMoves * 100.0);
            printf("Total energy is: %lf, error: %.15lf\n", state.energy, state.error);

            //1. Lista/vektor med olika input som de olika samplingsmetoderna behöver
            //2. sampler kan på något sätt efterfråga input, text genom att sätta en variabel
            //   Sen kan simulator ha en map och leta på den variabeln
            //sampler.sample(state);
            //for(auto s : sampler){
            //    s.sample(??????);
            //    s.sample(s.arguments);
            //}
            
        }
    }
};

#ifndef PY11
int main(){

    //trans.operator()<decltype(ps[1])>(ps[0]);

    Simulator* sim = new Simulator(2.0, 298.0);

    sim->state.set_geometry(0);
    sim->state.set_energy(0);
    /*std::vector< double > b;
    std::vector< double > q;
    b.push_back(0.0);
    b.push_back(0.0);
    q.push_back(1.0);
    q.push_back(-1.0);
    std::vector< std::vector<double> > pos;
    pos.emplace_back();
    pos.back() = {1, 2, 3};
    pos.emplace_back();
    pos.back() = {2, 2, 3};
    sim->state.particles.load(pos, q, b);*/
    sim->state.particles.create(200, 200);
    sim->state.equilibrate();
    sim->state.add_images();
    sim->state.finalize();
    sim->run(100, 1000);
    sim->state.particles.to_xyz("hej.xyz");
    //std::function<void(std::vector<int>)> move_callback = [state](std::vector<int> indices) { state.move_callback(indices); }

    return 0;
}
#endif

#ifdef PY11
PYBIND11_MODULE(mormon, m) {
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<int, int, double, double>())
        .def("run", &Simulator::run)
        .def_readwrite("state", &Simulator::state);

    py::class_<State>(m, "State")
        .def("load_state", &State::load_state)
        .def("set_geometry", &State::set_geometry)
        .def("set_energy", &State::set_energy)
        .def("load_particles", &State::load_particles)
        .def("equilibrate", &State::equilibrate)
        .def("finalize", &State::finalize);
}
#endif