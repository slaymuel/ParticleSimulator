#ifdef TRACK_MEMORY
    #include "memory_tracker.h"
#endif

#ifdef PY11
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    namespace py = pybind11;
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif

// Pre-compiled header
#include "pch.h"
#include "random/random.h"
#include "particles.h"
#include "state.h"
#include "move.h"
#include "sampler.h"
#include "comparators.h"
#include "io.h"
// timer defines the TIMEIT macro
#include "timer.h"

namespace Simulator{

class Simulator{
    
    private:
    //Name will be used for output files
    std::string name;
    //The relative probability weights of the moves
    std::vector<double> mWeights;
    //The list of moves
    std::vector< std::unique_ptr<Moves::Move> > moves;
    //List of samplers
    std::vector< std::unique_ptr<Samplers::SamplerBase> > sampler;

    public:
    //The current state of the system, i.e particles, geometry etc
    State state;

    Simulator(double Dielec, double T, std::string _name);
    // Reset/set the temperature, as in parallel tempering
    void set_temperature(double T);
    void add_move(Moves::MoveTypes move_type, std::vector<double> args);
    void add_sampler(Samplers::SamplerTypes type, std::vector<double> args);
    // Finalize the settings, called by user when user is done adding moves etc.
    void finalize();
    // Run the simulation
    void run(unsigned int macroSteps, unsigned int microSteps, unsigned int eqSteps);
};

// Python bindings using Pybind 11
#include "pybind.h"

} // end of namespace Simulator