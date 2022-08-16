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
    std::vector< std::unique_ptr<Move> > moves;
    std::vector< std::unique_ptr<Samplers::SamplerBase> > sampler;

    public:
    //The current state of the system
    State state;

    Simulator(double Dielec, double T, std::string _name);
    void set_temperature(double T);
    void set_cp(double cp);
    void add_move(MoveTypes move_type, std::vector<double> args);
    void add_sampler(int i, int interval, double ds = 0.05);
    //void add_sampler(Samplers::SamplerTypes type, std::vector<double> args);
    std::unique_ptr<Samplers::SamplerBase> _add_sampler(int i, int interval, double ds = 0.05);
    void finalize();
    void run(unsigned int macroSteps, unsigned int microSteps, unsigned int eqSteps);
};

#include "pybind.h"

}