#pragma once

namespace Simulator{

//Contains the global physical constants
namespace constants{
    constexpr double EC = 1.602176620898E-19;
    constexpr double VP = 8.854187817E-12;
    constexpr double PI = 3.14159265358979323846;
    constexpr double KB = 1.3806485279E-23;
    constexpr double C = EC * EC / (4.0 * PI * VP * 1e-10 * KB);
    extern double T;   //Temperature
    extern double cp;
    extern double D;   //Dielectric constant
    extern double lB;

    static void set(double _T, double _D){
        T = _T;
        D = _D;
        lB = C * (1.0 / (D * T));
    }
};

}