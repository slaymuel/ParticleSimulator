#pragma once

namespace Simulator{

namespace constants{
    static constexpr double EC = 1.602176620898E-19;
    static constexpr double VP = 8.854187817E-12;
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double KB = 1.3806485279E-23;
    static constexpr double C = EC * EC / (4.0 * PI * VP * 1e-10 * KB);
    static double T;   //Temperature
    static double cp;
    static double D;   //Dielectric constant
    static double lB;

    static void set(double _T, double _D){
        T = _T;
        D = _D;
        lB = C * (1.0 / (D * T));
    }
}

}