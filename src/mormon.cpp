#ifdef TRACK_MEMORY
    #include "memory_tracker.h"
#endif

#ifdef PY11
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    namespace py = pybind11;
#endif

#ifdef TIMERS
    #define TIMEIT Timer timer(__FUNCTION__);
#else
    #define TIMEIT;
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "pch.h"

#include "random.h"
#include "particles.h"
#include "state.h"
//#include <source_location>

#include "move.h"
#include "sampler.h"
#include "comparators.h"
#include "io.h"
#include "timer.h"




class Simulator{
    
    private:
    //Name will be used for output files
    std::string name;
    std::vector<double> mWeights;
    std::vector<double>::iterator wIt;
    std::vector<Move*> moves;
    std::vector<Samplers::SamplerBase*> sampler;


    public:
    State state;
    Simulator(double Dielec, double T, std::string name){
        constants::set(T, Dielec);

        this->name = name;

        #ifdef _OPENMP
            printf("\nOpenMP is ENABLED with %i threads.\n\n", omp_get_num_procs());
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
    
    void set_temperature(double T){
        constants::T = T; 
        constants::lB = constants::C * (1.0 / (constants::D * T));
    }
    void set_cp(double cp){
        constants::cp = cp;
    }


    void add_move(int i, double i1, double i2, double i3 = 0.0, double i4 = 0.0){
    //dp = i1, p = i2, cp = i3, d = i4
        printf("\nAdding move:\n");
        switch(i){
            case 0:
                moves.push_back(new Translate(i1, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 1:
                moves.push_back(new GrandCanonicalSingle<true>(i3, i4, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 2:
                moves.push_back(new GrandCanonicalSingle<false>(i3, i4, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 3:
                moves.push_back(new Rotate(i1, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 4:
                moves.push_back(new Swap(i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 5:
                moves.push_back(new SingleSwap(i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 6:
                moves.push_back(new VolumeMove(i1, i3, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 7:
                moves.push_back(new ChargeTrans(i1, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 8:
                moves.push_back(new ChargeTransRand(i1, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 9:
                moves.push_back(new Cluster(i1, i2, i3, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 10:
                moves.push_back(new WidomInsertion(i1, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 11:
                moves.push_back(new WidomDeletion(i1, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 12:
                moves.push_back(new GrandCanonical<true>(i3, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 13:
                moves.push_back(new GrandCanonical<false>(i3, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            case 14:
                moves.push_back(new ChargeTranslate(i1, i2, &state, std::bind(&State::move_callback, &state, std::placeholders::_1)));
                break;
            default:
                printf("Could not find move %i\n", i);
                break;
        }
    }

    void add_sampler(int i, int interval, double ds = 0.05){
        switch(i){
            case 0:
                printf("\nAdding z density sampler\n");
                sampler.push_back(new Samplers::Density(2, this->state.geo->_d[2], ds, 
                                              this->state.geo->d[0], this->state.geo->d[1], interval, this->name));
                break;
            case 1:
                printf("Adding Widom HS-CP sampler\n");
                sampler.push_back(new Samplers::WidomHS(interval, this->name));
                break;

            case 2:
                printf("Adding energy sampler\n");
                sampler.push_back(new Samplers::Energy(interval, this->name));
                break;

            case 3:
                printf("Adding charge distribution sampler\n");
                sampler.push_back(new Samplers::QDist(4, ds, interval, this->name));
                break;
            case 4:
                printf("Adding XDR trajectory sampler\n");
                sampler.push_back(new Samplers::XDR(interval, this->name));
                break;
            case 5:
                printf("Adding number of ions sampler\n");
                sampler.push_back(new Samplers::NumIons(interval, this->name));
                break;
            case 6:
                printf("\nAdding x density sampler\n");
                sampler.push_back(new Samplers::Density(0, this->state.geo->_d[0], ds, 
                                              this->state.geo->d[1], this->state.geo->d[2], interval, this->name));
                break;
            case 7:
                printf("\nAdding y density sampler\n");
                sampler.push_back(new Samplers::Density(1, this->state.geo->_d[1], ds, 
                                              this->state.geo->d[0], this->state.geo->d[2], interval, this->name));
                break;
            case 8:
                printf("\nAdding virial pressure sampler\n");
                sampler.push_back(new Samplers::Pressure(interval, this->state.geo->volume, this->state.geo->dh[2], this->name));
                break;
            case 9:
                printf("\nAdding pressureV sampler\n");
                sampler.push_back(new Samplers::PressureV(interval, ds, this->state.geo->_d[0], this->state.geo->_d[1], this->state.geo->_d[2], this->name));
                break;
            case 10:
                printf("\nAdding ForcePressure sampler\n");
                sampler.push_back(new Samplers::ForcePressure(interval, this->state.geo->volume, this->state.geo->dh[2], this->name));
                break;
            case 11:
                printf("\nAdding Force sampler\n");
                sampler.push_back(new Samplers::Force(interval, this->name));
                break;
            case 12:
                printf("\nAdding Cliff pressure sampler\n");
                sampler.push_back(new Samplers::CliffPressure(interval, ds, this->state.geo->_d[0], this->state.geo->_d[1], this->name));
                break;
            case 13:
                printf("\nAdding Modified Widom sampler\n");
                sampler.push_back(new Samplers::ModifiedWidom(interval, this->name));
                break;
            case 14:
                printf("\nAdding Modified Widom Coulomb sampler\n");
                sampler.push_back(new Samplers::ModifiedWidomCoulomb(interval, this->name));
                break;
            default:
                break;
        }
    }

    void finalize(){
        std::for_each( this->moves.begin(), this->moves.end(), [&](Move* m){ this->mWeights.push_back(m->weight); } );
        std::sort(this->moves.begin(), this->moves.end(), comparators::mLess);
        std::sort(this->mWeights.begin(), this->mWeights.end());

        for(unsigned int i = 1; i < this->mWeights.size(); i++){
            this->mWeights[i] += this->mWeights[i - 1];
        }

        //Make sure move list is not corrupted
        assert(this->mWeights.back() == 1.0);

        this->state.finalize(this->name);

        //Save starting configuration for XTC trajectory
        IO::to_gro(this->name, state.particles, state.geo->d);
    }

    void run(unsigned int macroSteps, unsigned int microSteps, unsigned int eqSteps){

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


        printf("Bjerrum length is: %.15lf\n", constants::lB);
        std::cout << "Running simulation at: " << constants::T << "K, "
                                                                << " with: " 
                                                                << state.particles.particles.size() << " particles ( "
                                                                << state.particles.cTot << " cations, " 
                                                                << state.particles.aTot <<" anions)" 
                                                                << std::endl;

        for(unsigned int macro = 0; macro < macroSteps; macro++){
            TIMEIT;
            for(unsigned int micro = 0; micro <= microSteps; micro++){
                wIt = std::lower_bound(mWeights.begin(), mWeights.end(), Random::get_random());
                (*moves[wIt - mWeights.begin()])();

                if(moves[wIt - mWeights.begin()]->accept( state.get_energy_change() )){
                    state.save();
                }
                else{
                    state.revert();
                }

                    //should also be able to
                    //state.get_energy(subset_of_particles);

                if(macro >= eqSteps){
                    for(auto s : sampler){
                        if(micro % s->interval == 0){
                            s->sample(state);  
                        }
                    }
                }
            }

            /*                                "HALF TIME"                                  */
            //Print progress
            std::cout << "\nIteration (macrostep): " << macro << std::endl;

            printf("Acceptance ratios: \n");
            for(auto move : moves){
                std::cout << move->dump() << std::endl;
            }
            
            printf("Total energy is: %lf, energy drift: %.15lf\n", state.cummulativeEnergy, state.error);
            printf("Cations: %i Anions: %i Tot: %i\n", state.particles.cTot, state.particles.aTot, state.particles.tot);
            printf("Box: %lf (%.15lf * %lf * %lf (%lf))\n", state.geo->volume, state.geo->_d[0], state.geo->_d[1], state.geo->_d[2], state.geo->d[2]);

            #ifdef TRACK_MEMORY
                std::cout << "Total allocated memory: " << allocationData.GetCurrentUsage() << std::endl << std::endl;
            #endif
            
            //Check energy drift etc
            state.control();
            state.advance();
            //1. Lista/vektor med olika input som de olika samplingsmetoderna behöver
            //2. sampler kan på något sätt efterfråga input, text genom att sätta en variabel
            //   Sen kan simulator ha en map och leta på den variabeln

            //for(auto s : sampler){
            //    s.sample(??????);
            //    s.sample(s.arguments);
            //}
            for(auto s : sampler){
                s->save();
            }
        }
        /*printf("Saving analysis data...\n");
        for(auto s : sampler){
            s->save(this->name);
        }*/

        for(auto s : sampler){
            s->close();
        }

        IO::to_xyz(this->name, state.particles, state.geo->d);
        IO::to_cpt(this->name, state.particles, state.geo->d);

        //this->state.close();
        printf("Energy of last frame: %.15lf\n", this->state.cummulativeEnergy);
        printf("Simulation Done!\n\n");
    }
};



#ifndef PY11
int main(){

    //trans.operator()<decltype(ps[1])>(ps[0]);
    std::string infile = "fgarpm_bulk.cp";
    std::string outfile = "fgarpm_bulk";
    
    double cutoff = 30.0;
    Simulator* sim = new Simulator(2.0, 900.0, outfile);

    sim->state.set_geometry(0, std::vector<double>{60.0, 60.0, 60.0});
    sim->state.set_energy(1, std::vector<double>{cutoff, 7, 7, 7, constants::PI / cutoff, 1, false});
    //sim->state.set_energy(9, {0.00425});

    //sim->state.particles.create(378, 378, 1.0, -1.0 , 2.5, 2.5, 2.5, 2.5, 0.0, 0.0);
    //sim->state.load_cp(infile);
    //sim->state.particles.set_models(std::vector<double>{1.0, -1.0}, std::vector<double>{2.5, 2.5}, std::vector<double>{2.5, 2.5},
    //                                std::vector<double>{0.0, 0.0}, std::vector<std::string>{"Na", "Cl"});
    sim->add_move(0, 0.12, 0.49);
    sim->add_move(0, 30.0, 0.01);
    sim->add_move(1, 0.0, 0.25, -16.0, 0.0);
    sim->add_move(2, 0.0, 0.25, -16.0, 0.0);
    //sim->add_move(6, 0.00025, 0.005, 101325.0);
    //sim->add_move(8, 0.2, 0.3);


    //sim->add_sampler(0, 100);
    //sim->add_sampler(2, 100);
    //sim->add_sampler(3, 100);
    //sim->add_sampler(5, 100);

    sim->state.equilibrate(10.0);

    //After equilibrate
    sim->finalize();
    sim->run(500, 10000, 0);
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
        .def("add_sampler", &Simulator::add_sampler, py::arg("i"), py::arg("interval"), py::arg("ds") = 0.05)
        .def("set_temperature", &Simulator::set_temperature)
        .def("set_cp", &Simulator::set_cp)
        .def("finalize", &Simulator::finalize)
        .def_readwrite("state", &Simulator::state);


    py::class_<State>(m, "State")
        .def("set_geometry", &State::set_geometry)
        .def("load_cp", &State::load_cp)
        .def("set_energy", &State::set_energy)
        .def("equilibrate", &State::equilibrate)
        .def("load_spline", &State::load_spline)
        .def("reset_energy", &State::reset_energy)
        .def_readwrite("particles", &State::particles)
        .def_readonly("energy", &State::energy)
        .def_readonly("cummulativeEnergy", &State::cummulativeEnergy);


    py::class_<Particles>(m, "Particles")
        .def_readonly("particles", &Particles::particles)
        .def_readonly("pModel", &Particles::pModel)
        .def_readonly("nModel", &Particles::nModel)
        //.def("load", &Particles::load)
        //.def("set_models", &Particles::set_models, py::arg("q"), py::arg("r"), py::arg("rf"), py::arg("b"), py::arg("names"))
        //.def("create", &Particles::create, py::arg("pNum"), py::arg("nNum"), py::arg("p"), py::arg("n"), py::arg("rfp") = 2.5, 
        //            py::arg("rfn") = 2.5, py::arg("rp") = 2.5, py::arg("rn") = 2.5, py::arg("bp") = 0.0, py::arg("bn") = 0.0);
        .def("create", &Particles::create, py::arg("pNum"), py::arg("nNum"), py::arg("params"));

    py::class_<Particle>(m, "Particle")
        .def_readonly("com", &Particle::com)
        .def_readonly("pos", &Particle::pos)
        .def_readonly("qDisp", &Particle::qDisp)
        .def_readonly("q", &Particle::q)
        .def_readonly("b", &Particle::b)
        .def_readonly("r", &Particle::r)
        .def_readonly("rf", &Particle::rf);
}
#endif