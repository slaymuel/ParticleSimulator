#pragma once

#include "aux_math.h"
#include "geometry.h"
#include "energy.h"
#include "potentials.h"
#include "Spline.h"

namespace Simulator{

// Forward declarations of friends
namespace Samplers{
    class PressureV;
    class CliffPressure;
}

class State{

    friend class Samplers::PressureV;
    friend class Samplers::CliffPressure;

    private:
    SplineData spline;
    std::shared_ptr<State> _old;
    //IO io;
    
    public:

    int step = 0;
    double energy = 0.0, cummulativeEnergy = 0.0, dE = 0.0, error = 0.0;
    Particles particles;
    std::vector< unsigned int > movedParticles;    //Particles that has moved from previous state
    std::shared_ptr<Geometry> geo;
    std::vector< std::shared_ptr<EnergyBase> > energyFunc;

    ~State(){}

    void advance();
    void control();
    void finalize(std::string name);
    void reset_energy();
    void save();
    void revert();
    //Get energy different between *this and old state
    double get_energy_change();
    //Called after move - set movedParticles
    void move_callback(std::vector< unsigned int > ps);
    void equilibrate(double step);
    bool overlap(std::size_t i) const;
    int get_overlaps() const;
    void set_geometry(int type, std::vector<double> args);
    void set_energy(int type, std::vector<double> args = std::vector<double>());
    void load_spline(std::vector<double> aKnots, std::vector<double> bKnots, std::vector<double >controlPoints);
    void load_cp(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, 
                 std::vector<double> charges, std::vector<double> r, std::vector<double> rf, 
                 std::vector<double> b, std::vector<double> b_min, std::vector<double> b_max, 
                 std::vector<std::string> names);
};

} // end of namespace Simulator