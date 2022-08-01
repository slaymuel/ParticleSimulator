#ifdef TRACK_MEMORY
    #include "memory_tracker.h"
#endif

#ifdef PY11
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    namespace py = pybind11;
#endif

#ifdef _TIMERS_
    #define TIMEIT Timer timer(__FUNCTION__);
#else
    #define TIMEIT;
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "pch.h"
#include "random/random.h"
#include "particles.h"
#include "state.h"
//#include <source_location>

#include "move.h"
#include "sampler.h"
#include "comparators.h"
#include "io.h"
#include "timer.h"

namespace Simulator{

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

    Simulator(double Dielec, double T, std::string _name) : name(_name){
        constants::set(T, Dielec);
        Move::state = &state;
    
        #ifdef _OPENMP
            Logger::Log("OpenMP is ENABLED with ", omp_get_num_procs()," threads.");
        #else
            Logger::Log("OpenMP is DISABLED");
        #endif

        #ifdef _DEBUG_
            Logger::Log<Logger::LogLevel::DEBUG>("Debug mode is ENABLED");
        #else
            Logger::Log("Debug mode is DISABLED");
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

    void add_move(MoveTypes moveType, std::vector<double> args){
        moves.push_back(Move::create_move(moveType, args));
    }

    void add_move(int i, std::vector<double> args){
    //dp = i1, p = i2, cp = i3, d = i4
        printf("\n");
        Logger::Log("Adding move:");
        switch(i){
            case 0:
                assert(args.size() == 2);
                                            //step, w
                moves.push_back(new Translate(args[1], args[0]));
                break;
            case 1:
                assert(args.size() == 3);
                                                                //cp, d, w
                moves.push_back(new GrandCanonicalSingle<true>(args[1], args[2], args[0]));
                break;
            case 2:
                assert(args.size() == 3);
                moves.push_back(new GrandCanonicalSingle<false>(args[1], args[2], args[0]));
                break;
            case 3:
                assert(args.size() == 2);
                                            //step, w
                moves.push_back(new Rotate(args[0], args[1]));
                break;
            case 4:
                assert(args.size() == 1);
                moves.push_back(new Swap(args[0]));
                break;
            case 5:
                assert(args.size() == 1);
                moves.push_back(new SingleSwap(args[0]));
                break;
            case 6:
                assert(args.size() == 3);
                                            //step, press, v, w
                moves.push_back(new VolumeMove(args[0], args[2], args[1]));
                break;
            case 7:
                assert(args.size() == 2);
                moves.push_back(new ChargeTrans(args[0], args[1]));
                break;
            case 8:
                assert(args.size() == 2);
                moves.push_back(new ChargeTransRand(args[0], args[1]));
                break;
            case 9:
                assert(args.size() == 3);
                moves.push_back(new Cluster(args[0], args[1], args[2]));
                break;
            case 10:
                assert(args.size() == 1);
                moves.push_back(new WidomInsertion(args[0]));
                break;
            case 11:
                assert(args.size() == 1);
                moves.push_back(new WidomDeletion(args[0]));
                break;
            case 12:
                assert(args.size() == 2);
                                                        //chempot, w
                moves.push_back(new GrandCanonical<true>(args[1], args[0]));
                break;
            case 13:
                assert(args.size() == 2);
                moves.push_back(new GrandCanonical<false>(args[1], args[0]));
                break;
            case 14:
                assert(args.size() == 2);
                moves.push_back(new ChargeTranslate(args[0], args[1]));
                break;
            default:
                Logger::Log("Could not find move ", i, "\n");
                break;
        }
    }

    void add_sampler(int i, int interval, double ds = 0.05){
        switch(i){
            case 0:
                Logger::Log("\nAdding z density sampler\n");
                sampler.push_back(new Samplers::Density(2, this->state.geo->_d[2], ds, 
                                              this->state.geo->d[0], this->state.geo->d[1], interval, this->name));
                break;
            case 1:
                Logger::Log("Adding Widom HS-CP sampler\n");
                sampler.push_back(new Samplers::WidomHS(interval, this->name));
                break;

            case 2:
                Logger::Log("Adding energy sampler\n");
                sampler.push_back(new Samplers::Energy(interval, this->name));
                break;

            case 3:
                Logger::Log("Adding charge distribution sampler\n");
                sampler.push_back(new Samplers::QDist(4, ds, interval, this->name));
                break;
            case 4:
                Logger::Log("Adding XDR trajectory sampler\n");
                sampler.push_back(new Samplers::XDR(interval, this->name));
                break;
            case 5:
                Logger::Log("Adding number of ions sampler\n");
                sampler.push_back(new Samplers::NumIons(interval, this->name));
                break;
            case 6:
                Logger::Log("\nAdding x density sampler\n");
                sampler.push_back(new Samplers::Density(0, this->state.geo->_d[0], ds, 
                                              this->state.geo->d[1], this->state.geo->d[2], interval, this->name));
                break;
            case 7:
                Logger::Log("\nAdding y density sampler\n");
                sampler.push_back(new Samplers::Density(1, this->state.geo->_d[1], ds, 
                                              this->state.geo->d[0], this->state.geo->d[2], interval, this->name));
                break;
            case 8:
                Logger::Log("\nAdding virial pressure sampler\n");
                sampler.push_back(new Samplers::Pressure(interval, this->state.geo->volume, this->state.geo->dh[2], this->name));
                break;
            case 9:
                Logger::Log("\nAdding pressureV sampler\n");
                sampler.push_back(new Samplers::PressureV(interval, ds, this->state.geo->_d[0], this->state.geo->_d[1], this->state.geo->_d[2], this->name));
                break;
            case 10:
                Logger::Log("\nAdding ForcePressure sampler\n");
                sampler.push_back(new Samplers::ForcePressure(interval, this->state.geo->volume, this->state.geo->dh[2], this->name));
                break;
            case 11:
                Logger::Log("\nAdding Force sampler\n");
                sampler.push_back(new Samplers::Force(interval, this->name));
                break;
            case 12:
                Logger::Log("\nAdding Cliff pressure sampler\n");
                sampler.push_back(new Samplers::CliffPressure(interval, ds, this->state.geo->_d[0], this->state.geo->_d[1], this->name));
                break;
            case 13:
                Logger::Log("\nAdding Modified Widom sampler\n");
                sampler.push_back(new Samplers::ModifiedWidom(interval, this->name));
                break;
            case 14:
                Logger::Log("\nAdding Modified Widom Coulomb sampler\n");
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
        Logger::Log("Bjerrum length is: ", constants::lB);
        Logger::Log("Running Simulation at: ", constants::T, "K", " with ", state.particles.particles.size(), " particles (", state.particles.cTot, " cations and ", state.particles.aTot, " anions)");

        for(unsigned int macro = 0; macro < macroSteps; macro++){
            #ifdef _TIMERS_
                TIMEIT;
            #endif
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
            printf("\n");
            Logger::Log("Iteration (macrostep): ", macro);

            Logger::Log("Acceptances ratios");
            for(auto move : moves){
                Logger::Log(move->dump());
            }
            
            Logger::Log("Total energy is: ", state.cummulativeEnergy, ", energy drift: ", state.error);
            Logger::Log("Cations: ", state.particles.cTot, " Anions: ", state.particles.aTot, " Tot: ", state.particles.tot);
            Logger::Log("Box: ", state.geo->volume, " (", state.geo->_d[0] ," * " , state.geo->_d[1], " * ", state.geo->_d[2], " (", state.geo->d[2], ") )");

            #ifdef TRACK_MEMORY
                Logger::Log("Total allocated memory: ", allocationData.GetCurrentUsage());
            #endif
            
            //Check energy drift etc
            state.control();
            state.advance();

            for(auto s : sampler){
                s->save();
            }
        }

        for(auto s : sampler){
            s->close();
        }

        IO::to_xyz(this->name, state.particles, state.geo->d);
        IO::to_cpt(this->name, state.particles, state.geo->d);

        //this->state.close();
        Logger::Log("Energy of last frame: ", this->state.cummulativeEnergy);
        Logger::Log("Simulation Done!");
        printf("\n\n");
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
PYBIND11_MODULE(particlesimulator, m) {
    py::enum_<MoveTypes>(m, "MoveTypes")
        .value("TRANSLATE", MoveTypes::TRANSLATE)
        .value("GCSINGLEADD", MoveTypes::GCSINGLEADD)
        .value("GCSINGLEREMOVE", MoveTypes::GCSINGLEREMOVE)
        .value("ROTATE", MoveTypes::ROTATE)
        .value("SWAP", MoveTypes::SWAP)
        .value("SINGLESWAP", MoveTypes::SINGLESWAP)
        .value("VOLUMEMOVE", MoveTypes::VOLUMEMOVE)
        .value("CHARGETRANS", MoveTypes::CHARGETRANS)
        .value("CHARGETRANSRAND", MoveTypes::CHARGETRANSRAND)
        .value("CLUSTER", MoveTypes::CLUSTER)
        .value("WIDOMINSERTION", MoveTypes::WIDOMINSERTION)
        .value("WIDOMDELETION", MoveTypes::WIDOMDELETION)
        .value("GCADD", MoveTypes::GCADD)
        .value("GCREMOVE", MoveTypes::GCREMOVE)
        .value("CHARGETRANSLATE", MoveTypes::CHARGETRANSLATE);

    

    py::class_<Simulator>(m, "Simulator")
        .def(py::init<double, double, std::string>())
        .def("run", &Simulator::run)
        //.def("add_move", &Simulator::add_move)
        .def("add_move", static_cast<void (Simulator::*)(MoveTypes, std::vector<double>)>(&Simulator::add_move))
        .def("add_move",  static_cast<void (Simulator::*)(int, std::vector<double>)>(&Simulator::add_move))
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

}