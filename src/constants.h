#pragma once

namespace Simulator{

//Contains the global physical constants
namespace constants{
    // Elementary charge
    inline constexpr double EC = 1.602176620898E-19;
    // Vacuum permittivity
    inline constexpr double VP = 8.854187817E-12;
    // Pi
    inline constexpr double PI = 3.14159265358979323846;
    // Boltzmann constant
    inline constexpr double KB = 1.3806485279E-23;
    inline constexpr double C = EC * EC / (4.0 * PI * VP * 1e-10 * KB);
    // Temperature
    extern double T;
    
    extern double D;
    // Bjerrum length
    extern double lB;

    inline void set(double _T, double _D){
        T = _T;
        D = _D;
        lB = C * (1.0 / (D * T));
    }
}

} // end of namespace Simulator