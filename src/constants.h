#pragma once

namespace Simulator{

//Contains the global physical constants
namespace constants{
    inline constexpr double EC = 1.602176620898E-19;
    inline constexpr double VP = 8.854187817E-12;
    inline constexpr double PI = 3.14159265358979323846;
    inline constexpr double KB = 1.3806485279E-23;
    inline constexpr double C = EC * EC / (4.0 * PI * VP * 1e-10 * KB);
    extern double T;   //Temperature
    extern double cp;
    extern double D;   //Dielectric constant
    extern double lB;

    inline void set(double _T, double _D){
        T = _T;
        D = _D;
        lB = C * (1.0 / (D * T));
    }
}

} // end of namespace Simulator