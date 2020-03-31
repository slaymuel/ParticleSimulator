#pragma once

namespace constants{
    const double EC = 1.602176620898E-19;
    const double VP = 8.854187817E-12;
    const double PI = 3.14159265358979323846;
    const double KB = 1.3806485279E-23;
    const double C = EC * EC / (4.0 * PI * VP * 1e-10 * KB);
    double T;   //Temperature
    double cp;
    double D;   //Dielectric constant
    double lB;
}

