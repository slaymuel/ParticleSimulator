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
    // The old state
    std::shared_ptr<State> _old;
    
    public:
    // The current step
    int step = 0;
    // Energy variables
    double energy = 0.0, cummulativeEnergy = 0.0, dE = 0.0, error = 0.0;
    // The particles object which holds all particles that are being simulated
    Particles particles;
    // Populated by move_callback, holds moved particles
    std::vector< unsigned int > movedParticles;
    // The geometry used in the simulation
    std::shared_ptr<Geometry> geo;
    // Holds all the energy functors
    std::vector< std::shared_ptr<EnergyBase> > energyFunc;

    ~State(){}

    // Advance the simulation one step
    void advance();
    // Control mechanisms for debugging
    void control();
    // Called by the user when the settings are finished
    void finalize(std::string name);
    // Recalculate the energy and reset the cummulative energy
    void reset_energy();
    // Save the state, called when a move is accepted
    void save();
    // Revert the state, i.e state = old, called when a move is rejected
    void revert();
    //Get energy different between *this and old state
    double get_energy_change();
    //Called after move - set movedParticles
    void move_callback(std::vector< unsigned int > ps);
    // Equilibrate the system, i.e remove overlaps
    void equilibrate(double step);
    //Checks if the particle at index i overlaps with any other particle
    bool overlap(std::size_t i) const;
    // Calculate the number of overlaps
    int get_overlaps() const;
    // Set the simulation geometry
    void set_geometry(int type, std::vector<double> args);
    // Add a potential
    void set_energy(int type, std::vector<double> args = std::vector<double>());
    // Load a potential spline from file
    void load_spline(std::vector<double> aKnots, std::vector<double> bKnots, std::vector<double >controlPoints);
    // Loads a data from a checkpoint
    void load_cp(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, 
                 std::vector<double> charges, std::vector<double> r, std::vector<double> rf, 
                 std::vector<double> b, std::vector<double> b_min, std::vector<double> b_max, 
                 std::vector<std::string> names);
};

} // end of namespace Simulator