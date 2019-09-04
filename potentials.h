#include "particle.h"

class Coulomb{
    public:

    inline double operator()(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, double dist){
        return p1->q * p2->q / dist;
    }
};

//In GC ewald should only return reciprocal part in previous state