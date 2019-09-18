#include <random>
#include "random.h"
#include <assert.h>
#include "constants.h"
#include <iostream>
#include <omp.h>
#include "aux_math.h"
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
    std::vector<double>::iterator wIt;
    std::vector<Move*> moves;
    std::vector<Sampler*> sampler;
    /* State callback after move */
    std::function< void(std::vector< unsigned int >) > move_callback 
                = std::bind(&State::move_callback, &state, std::placeholders::_1);
    std::string name;

    public:
    State state;
    Simulator(double Dielec, double T, std::string name){
        //Set some constants
        constants::D = Dielec;
        constants::T = T; 
        constants::lB = constants::C * (1.0 / (Dielec * T));
        this->name = name;
        Random::initialize();

        #ifdef _OPENMP
            printf("\nOpenMP is ENABLED with %i cores.\n\n", omp_get_num_procs());
        #else
            printf("\nOpenMP is DISABLED\n");
        #endif

        #ifdef DEBUG
        printf("Debug mode is ENABLED\n");
        #else
        printf("Debug mode is DISABLED\n");
        #endif
        printf("\n");
    }
    
    


    void add_move(int i, double dp, double p, double cp = 0.0, double d = 0.0){
        switch(i){
            case 0:
                moves.push_back(new Translate(dp, p));
                break;
            case 1:
                moves.push_back(new GrandCanonicalAdd(cp, d, &state, p));
                break;
            case 2:
                moves.push_back(new GrandCanonicalRemove(cp, d, &state, p));
                break;
            default:
                printf("Could not find move %i\n", i);
                break;
        }
    }

    void add_sampler(int i){
        switch(i){
            case 0:
                sampler.push_back(new Density(2, this->state.geo->d[2] / 2.0, 0.05, 
                                              this->state.geo->d[0], this->state.geo->d[1]));
        }
    }

    void run(int macroSteps, int microSteps){

        //////////////////// MOVE TO SOMEWHERE ELSE /////////////////////////////////////////////////////////
        std::for_each( moves.begin(), moves.end(), [&](Move* m){ mWeights.push_back(m->weight); } );

        std::sort(moves.begin(), moves.end(), comparators::mLess);
        std::sort(mWeights.begin(), mWeights.end());
        for(int i = 1; i < mWeights.size(); i++){
            mWeights[i] += mWeights[i - 1];
        }
        assert(mWeights.back() == 1.0);
        ////////////////////////////////////////////////////////////////////////////////////////////////////


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



        std::cout << "Running simulation at: " << constants::T << "K, lB: "<< constants::lB << " with: " << state.particles.particles.size() 
                                                                << " particles" << std::endl;

        for(int macro = 0; macro < macroSteps; macro++){
            auto start = std::chrono::steady_clock::now();
            for(int micro = 0; micro < microSteps; micro++){
                
                wIt = std::lower_bound(mWeights.begin(), mWeights.end(), Random::get_random());

                (*moves[wIt - mWeights.begin()])(state.particles.random(), move_callback);
                if(moves[wIt - mWeights.begin()]->accept( state.get_energy_change() )){
                    state.save();
                }
                else{
                    state.revert();
                }

                    //should also be able to
                    //state.get_energy(subset_of_particles);
                if(micro % 100 == 0 && micro > 0 && macro > 10){
                    for(auto s : sampler){
                        s->sample(state.particles);
                    }
                }
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
            std::cout << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << "s per macrostep\n\n";

            //1. Lista/vektor med olika input som de olika samplingsmetoderna behöver
            //2. sampler kan på något sätt efterfråga input, text genom att sätta en variabel
            //   Sen kan simulator ha en map och leta på den variabeln

            //for(auto s : sampler){
            //    s.sample(??????);
            //    s.sample(s.arguments);
            //}
        }
        printf("Saving analysis data...\n");
        for(auto s : sampler){
            s->save(this->name);
        }
        this->state.particles.to_xyz(this->name);
        printf("Simulation Done!\n\n");
    }
};



#ifndef PY11
int main(){

    //trans.operator()<decltype(ps[1])>(ps[0]);

    Simulator* sim = new Simulator(78.3, 298.0, "halfwald_test");
    sim->add_move(0, 5.0, 0.99);
    sim->add_move(1, 0.0, 0.005, -10.7, -2.0);
    sim->add_move(2, 0.0, 0.005, -10.7, -2.0);
    sim->state.set_geometry(2);
    sim->state.set_energy(2);
    sim->add_sampler(0);
    /*std::vector< double > b;
    std::vector< double > q;
    std::vector< std::string > n;
    n.push_back("Na");
    n.push_back("Cl");
    b.push_back(0.0);
    b.push_back(0.0);
    q.push_back(1.0);
    q.push_back(-1.0);
    std::vector< std::vector<double> > pos;
    pos.emplace_back();
    pos.back() = {0, 0, 4.5};
    pos.emplace_back();
    pos.back() = {0, 0, -4.5};
    sim->state.particles.load(pos, q, b, n);*/
                              // +    -
    sim->state.particles.create(373, 243, 2.0, -1.0);
    sim->state.equilibrate();
    //sim->state.add_images();
    sim->state.finalize();
    sim->run(100, 10000);
    //sim->state.particles.to_xyz("hej.xyz");
    //std::function<void(std::vector<int>)> move_callback = [state](std::vector<int> indices) { state.move_callback(indices); }

    return 0;
}
#endif

#ifdef PY11
PYBIND11_MODULE(mormon, m) {
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<double, double, std::string>())
        .def("run", &Simulator::run)
        .def("add_move", &Simulator::add_move, py::arg("i"), py::arg("dp"), py::arg("p"), py::arg("cp") = 0.0, py::arg("d") = 0.0)
        .def("add_sampler", &Simulator::add_sampler)
        .def_readwrite("state", &Simulator::state);

    py::class_<State>(m, "State")
        .def("load_state", &State::load_state)
        .def("set_geometry", &State::set_geometry)
        .def("set_energy", &State::set_energy)
        .def("equilibrate", &State::equilibrate)
        .def("finalize", &State::finalize)
        .def_readwrite("particles", &State::particles);

    py::class_<Particles>(m, "Particles")
        .def("load", &Particles::load)
        .def("create", &Particles::create, py::arg("pNum"), py::arg("nNum"), py::arg("p"), py::arg("n"), py::arg("rfp") = 2.5, py::arg("rfn") = 2.5);
}
#endif