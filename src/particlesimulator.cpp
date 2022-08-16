#include "particlesimulator.h"

namespace Simulator{

Simulator::Simulator(double Dielec, double T, std::string _name) : name(_name){
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

void Simulator::set_temperature(double T){
    constants::T = T; 
    constants::lB = constants::C * (1.0 / (constants::D * T));
}

void Simulator::set_cp(double cp){
    constants::cp = cp;
}

void Simulator::add_move(MoveTypes move_type, std::vector<double> args){
    moves.push_back(Move::createMove(move_type, args));
}

void Simulator::add_sampler(int i, int interval, double ds){
    sampler.push_back(_add_sampler(i, interval, ds));
}

/*void Simulator::add_sampler(Samplers::SamplerTypes type, std::vector<double> args){
    Samplers::createSampler(type, args, this->state.geo);
}*/

std::unique_ptr<Samplers::SamplerBase> Simulator::_add_sampler(int i, int interval, double ds){
    switch(i){
        case 0:
            Logger::Log("\nAdding z density sampler\n");
            return std::make_unique<Samplers::Density>(2, this->state.geo->_d[2], ds, 
                                            this->state.geo->d[0], this->state.geo->d[1], interval, this->name);
            break;
        case 1:
            Logger::Log("Adding Widom HS-CP sampler\n");
            return std::make_unique<Samplers::WidomHS>(interval, this->name);
            break;

        case 2:
            Logger::Log("Adding energy sampler\n");
            return std::make_unique<Samplers::Energy>(interval, this->name);
            break;

        case 3:
            Logger::Log("Adding charge distribution sampler\n");
            return std::make_unique<Samplers::QDist>(4, ds, interval, this->name);
            break;
        case 4:
            Logger::Log("Adding XDR trajectory sampler\n");
            return std::make_unique<Samplers::XDR>(interval, this->name);
            break;
        case 5:
            Logger::Log("Adding number of ions sampler\n");
            return std::make_unique<Samplers::NumIons>(interval, this->name);
            break;
        case 6:
            Logger::Log("\nAdding x density sampler\n");
            return std::make_unique<Samplers::Density>(0, this->state.geo->_d[0], ds, 
                                            this->state.geo->d[1], this->state.geo->d[2], interval, this->name);
            break;
        case 7:
            Logger::Log("\nAdding y density sampler\n");
            return std::make_unique<Samplers::Density>(1, this->state.geo->_d[1], ds, 
                                            this->state.geo->d[0], this->state.geo->d[2], interval, this->name);
            break;
        case 8:
            Logger::Log("\nAdding virial pressure sampler\n");
            return std::make_unique<Samplers::Pressure>(interval, this->state.geo->volume, this->state.geo->dh[2], this->name);
            break;
        case 9:
            Logger::Log("\nAdding pressureV sampler\n");
            return std::make_unique<Samplers::PressureV>(interval, ds, this->state.geo->_d[0], this->state.geo->_d[1], this->state.geo->_d[2], this->name);
            break;
        case 10:
            Logger::Log("\nAdding ForcePressure sampler\n");
            return std::make_unique<Samplers::ForcePressure>(interval, this->state.geo->volume, this->state.geo->dh[2], this->name);
            break;
        case 11:
            Logger::Log("\nAdding Force sampler\n");
            return std::make_unique<Samplers::Force>(interval, this->name);
            break;
        case 12:
            Logger::Log("\nAdding Cliff pressure sampler\n");
            return std::make_unique<Samplers::CliffPressure>(interval, ds, this->state.geo->_d[0], this->state.geo->_d[1], this->name);
            break;
        case 13:
            Logger::Log("\nAdding Modified Widom sampler\n");
            return std::make_unique<Samplers::ModifiedWidom>(interval, this->name);
            break;
        case 14:
            Logger::Log("\nAdding Modified Widom Coulomb sampler\n");
            return std::make_unique<Samplers::ModifiedWidomCoulomb>(interval, this->name);
            break;
        default:
            Logger::Log("Invalid sampler");
            break;
    }
}

void Simulator::finalize(){
    std::for_each( this->moves.begin(), this->moves.end(), 
                    [&](std::unique_ptr<Move>& m){ this->mWeights.push_back(m->weight); } );
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

void Simulator::run(unsigned int macroSteps, unsigned int microSteps, unsigned int eqSteps){
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

            if(macro >= eqSteps){
                for(auto& s : sampler){
                    if(micro % s->getInterval() == 0){
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
        for(const auto& move : moves){
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

        for(const auto& s : sampler){
            s->save();
        }
    }

    for(const auto& s : sampler){
        s->close();
    }

    IO::to_xyz(this->name, state.particles, state.geo->d);
    IO::to_cpt(this->name, state.particles, state.geo->d);

    //this->state.close();
    Logger::Log("Energy of last frame: ", this->state.cummulativeEnergy);
    Logger::Log("Simulation Done!");
    printf("\n\n");
}

#ifndef PY11
int main(){
    std::string infile = "fgarpm_bulk.cp";
    std::string outfile = "fgarpm_bulk";
    
    double cutoff = 30.0;
    Simulator sim(2.0, 900.0, outfile);

    sim.state.set_geometry(0, std::vector<double>{60.0, 60.0, 60.0});
    sim.state.set_energy(1, std::vector<double>{cutoff, 7, 7, 7, constants::PI / cutoff, 1, false});

    //sim->state.load_cp(infile);

    sim.add_move(0, 0.12, 0.49);
    sim.add_move(0, 30.0, 0.01);
    sim.add_move(1, 0.0, 0.25, -16.0, 0.0);
    sim.add_move(2, 0.0, 0.25, -16.0, 0.0);

    //sim->add_sampler(0, 100);

    sim.state.equilibrate(10.0);

    //After equilibrate
    sim.finalize();
    sim.run(500, 10000, 0);

    return 0;
}
#endif

} // end of namespace Simulator