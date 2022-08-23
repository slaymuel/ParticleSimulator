#include "particlesimulator.h"

namespace Simulator{

Simulator::Simulator(double Dielec, double T, std::string _name) : name(_name){
    constants::set(T, Dielec);
    Moves::Move::state = &state;

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

// Set the temperature
void Simulator::set_temperature(double T){
    constants::T = T; 
    constants::lB = constants::C * (1.0 / (constants::D * T));
}

void Simulator::add_move(Moves::MoveTypes move_type, std::vector<double> args){
    moves.push_back(Moves::moveFactory.createObject(move_type, this->name, args));
}

void Simulator::add_sampler(Samplers::SamplerTypes type, std::vector<double> args){
    //      0       1    2     3     4     5      6       7
    // { interval, ds, d[0], d[1], d[2], _d[0], _d[1], _d[2] }
    assert( (args.size() > 0) && (args.size() < 3) );

    // If user did not supply a ds argument, push back default value
    if(args.size() == 1)
        args.push_back(0.01);
    // Add values needed by the samplers
    args.insert(args.end(), this->state.geo->d.begin(), this->state.geo->d.end());
    args.insert(args.end(), this->state.geo->_d.begin(), this->state.geo->_d.end());
    sampler.push_back(Samplers::samplerFactory.createObject(type, this->name, args));
}

void Simulator::finalize(){
    // Create the vector with weights
    std::for_each( this->moves.begin(), this->moves.end(), 
                    [&](std::unique_ptr<Moves::Move>& m){ this->mWeights.push_back(m->weight); } );
    std::sort(this->moves.begin(), this->moves.end(), comparators::mLess);
    std::sort(this->mWeights.begin(), this->mWeights.end());

    // Calculate the accumulated weights
    for(unsigned int i = 1; i < this->mWeights.size(); i++)
        this->mWeights[i] += this->mWeights[i - 1];

    // Make sure move list is not corrupted
    assert(this->mWeights.back() == 1.0);

    // Sets up the old system and calculated the starting energy
    this->state.finalize(this->name);

    // Save starting configuration for XTC trajectory
    IO::instance().to_gro(this->name, state.particles, state.geo->d);
}

// Runner
void Simulator::run(unsigned int macroSteps, unsigned int microSteps, unsigned int eqSteps){
    Logger::Log("Bjerrum length is: ", constants::lB);
    Logger::Log("Running Simulation at: ", constants::T, "K", " with ", 
                state.particles.tot, " particles (", state.particles.cTot, 
                " cations and ", state.particles.aTot, " anions)");

    // Main program loop
    // For each macrostep
    for(unsigned int macro = 0; macro < macroSteps; macro++){
        #ifdef _TIMERS_
            TIMEIT;
        #endif
        
        // For each microstep
        for(unsigned int micro = 0; micro <= microSteps; micro++){
            // Get a random move using move probability distribution 
            auto wIt = std::lower_bound(mWeights.begin(), mWeights.end(), Random::get_random());
            // Perform the move
            (*moves[wIt - mWeights.begin()])();
            // If the move is accepted, save the state
            if(moves[wIt - mWeights.begin()]->accept( state.get_energy_change() ))
                state.save();
            // If the move is rejected, revert to the old state
            else
                state.revert();
            // Sample the state
            if(macro >= eqSteps){
                for(auto& s : sampler){
                    if(micro % s->getInterval() == 0)
                        s->sample(state);
                }
            }
            // Move to next micro step
            state.advanceMicro();
        } // end of inner microstepping loop

        // Print progress
        printf("\n");
        Logger::Log("Iteration (macrostep): ", macro);

        Logger::Log("Acceptances ratios");
        for(const auto& move : moves){
            Logger::Log(move->dump());
        }
        
        Logger::Log("Total energy is: ", state.cummulativeEnergy, ", energy drift: ", state.error);
        Logger::Log("Cations: ", state.particles.cTot, " Anions: ", state.particles.aTot, " Tot: ", 
                                                                                state.particles.tot);
        Logger::Log("Box: ", state.geo->volume, " (", state.geo->_d[0] ," * " , state.geo->_d[1], " * ", 
                                                        state.geo->_d[2], " (", state.geo->d[2], ") )");

        #ifdef TRACK_MEMORY
            Logger::Log("Total allocated memory: ", allocationData.GetCurrentUsage());
        #endif
        
        //Check energy drift etc
        state.control();
        
        // Save sampled data to file
        for(const auto& s : sampler)
            s->save();

        // Advance to the next macrostep
        state.advance();

    } // End of main program loop

    for(const auto& s : sampler)
        s->close();

    // Save final configuration and checkpoint file
    IO::instance().to_xyz(this->name, state.particles, state.geo->d);
    IO::instance().to_cpt(this->name, state.particles, state.geo->d);

    Logger::Log("Energy of last frame: ", this->state.cummulativeEnergy);
    Logger::Log("Simulation Done!");
    printf("\n\n");
}

// Deprecated, currently cannot use without Pybind11
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