#include <random>
#include "random.h"
#include <assert.h>
#include "constants.h"
#include <iostream>
#include "state.h"
#include "particle.h"
#include "move.h"
#include <functional>
#include <chrono>
#include "sampler.h"
#include <algorithm>
#include "comparators.h"
#ifdef PY11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
#endif



class Simulator{
    private:
    std::vector<Particle*> ps;
    std::vector<double> mWeights;
    public:
    State state;
    Simulator(double Dielec, double T){
        //Set some constants
        constants::D = Dielec;
        constants::T = T; 
        constants::lB = constants::C * (1.0 / (Dielec * T));
        Random::initialize();

        #ifdef _OPENMP
            printf("OpenMP is ENABLED\n");
        #else
            printf("OpenMP is DISABLED\n");
        #endif
    }

    //Energy<Coulomb> energy;
    //moves
    
    std::vector<Move*> moves;
    
    /* State callback after move */
    std::function< void(std::vector< int >) > move_callback 
                = std::bind(&State::move_callback, &state, std::placeholders::_1);


    void run(int macroSteps, int microSteps){
        printf("            +\n"                                            
               "           (|)\n"
               " _____.___.|_|.\n"                                      
               "|    / \\  |===|\n"                                           
               "|   /   \\ | o |               MORMONS\n"                                           
               "|__/__v__\\|, ,|    MOleculaR MOdelliNg Software\n"                                           
               "| | | | | || ||   -----------------------------\n"   
               "|/| . . . |','|\n"                                           
               "||| A A A | , |\n"                                           
               "||| M M M |   |\n"                                       
               "---------------\n\n");



        std::cout << "Running simulation at: " << constants::T << "K with: " << state.particles.particles.size() 
                                                                << " particles" << std::endl;


        //////////////////// MOVE TO SOMEWHERE ELSE /////////////////////////////////////////////////////////
        moves.push_back(new Translate<false>(1.0, 0.98));
        moves.push_back(new GrandCanonicalAdd<false>(-5.0, 0.0, &state, 0.01));
        moves.push_back(new GrandCanonicalRemove<false>(-5.0, 0.0, &state, 0.01));
        std::for_each( moves.begin(), moves.end(), [&](Move* m){ mWeights.push_back(m->weight); } );

        std::sort(moves.begin(), moves.end(), comparators::mLess);
        std::sort(mWeights.begin(), mWeights.end());
        for(int i = 1; i < mWeights.size(); i++){
            mWeights[i] += mWeights[i - 1];
        }
        assert(mWeights.back() == 1.0);
        std::vector<double>::iterator wIt;
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        Sampler* sampler = new Density(2, state.geo->d[2], 0.2);

        for(int macro = 0; macro < macroSteps; macro++){
            auto start = std::chrono::steady_clock::now();
            for(int micro = 0; micro < microSteps; micro++){
                
                wIt = std::lower_bound(mWeights.begin(), mWeights.end(), Random::get_random());

                //std::cout << "lower_bound at position " << (wIt - mWeights.begin()) << '\n';
                //Move should check if particle is part of molecule
                //(*moves[0])(state.particles.random(), move_callback); // Two virtual calls

                (*moves[wIt - mWeights.begin()])(state.particles.random(), move_callback);
                if(moves[wIt - mWeights.begin()]->accept( state.get_energy_change() )){
                    state.save();
                }
                else{
                    state.revert();
                }

                    //should also be able to
                    //state.get_energy(subset_of_particles);
            }
            /*                                HALF TIME                                  */
            //Check energy drift etc
            state.control();

            //Print progress
            std::cout << "\nIteration: " << macro * microSteps + microSteps << std::endl;

            printf("Acceptance ratios: ");
            for(auto move : moves){
                printf("%s %.1lf%% %i ", move->id.c_str(), (double)move->accepted / move->attempted * 100.0, move->attempted);
            }
            printf("\n");
            
            printf("Total energy is: %lf, energy drift: %.15lf\n", state.energy, state.error);
            printf("Cations: %i Anions: %i Tot: %i\n", state.particles.cTot, state.particles.aTot, state.particles.tot);
            auto end = std::chrono::steady_clock::now();
            std::cout << (double) std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / microSteps << "us\n\n";

            sampler->sample(state.particles);
            //1. Lista/vektor med olika input som de olika samplingsmetoderna behöver
            //2. sampler kan på något sätt efterfråga input, text genom att sätta en variabel
            //   Sen kan simulator ha en map och leta på den variabeln

            //for(auto s : sampler){
            //    s.sample(??????);
            //    s.sample(s.arguments);
            //}
        }
        printf("Saving analysis data...\n");
        sampler->save("z_dens.txt");
        printf("Simulation Done!\n\n");
    }
};

#ifndef PY11
int main(){

    //trans.operator()<decltype(ps[1])>(ps[0]);

    Simulator* sim = new Simulator(78.0, 298.0);

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
    sim->state.particles.create(100, 100);
    sim->state.equilibrate();
    //sim->state.add_images();
    sim->state.finalize();
    sim->run(1000, 1000);
    sim->state.particles.to_xyz("hej.xyz");
    //std::function<void(std::vector<int>)> move_callback = [state](std::vector<int> indices) { state.move_callback(indices); }

    return 0;
}
#endif

#ifdef PY11
PYBIND11_MODULE(mormon, m) {
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<double, double>())
        .def("run", &Simulator::run)
        .def_readwrite("state", &Simulator::state);

    py::class_<State>(m, "State")
        .def("load_state", &State::load_state)
        .def("set_geometry", &State::set_geometry)
        .def("set_energy", &State::set_energy)
        .def("equilibrate", &State::equilibrate)
        .def("finalize", &State::finalize)
        .def_readwrite("particles", &State::particles);

    py::class_<Particles>(m, "Particles")
        .def("create", &Particles::create);
}
#endif