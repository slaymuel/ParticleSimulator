#pragma once

#include "state.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"
#include "genericfactory.h"

namespace Simulator{

namespace Samplers{

class SamplerBase;

enum class SamplerTypes{
    DENSITY_X,
    DENSITY_Y,
    DENSITY_Z,
    ENERGY,
    WIDOMHS,
    QDIST,
    XDR,
    NUMIONS,
    PRESSURE,
    PRESSUREV,
    FORCEPRESSURE,
    CLUSTER,
    FORCE,
    CLIFFPRESSURE,
    MODIFIEDWIDOM,
    MODIFIEDWIDOMCOULOMB
};

//      0       1    2     3     4     5      6       7
// { interval, ds, d[0], d[1], d[2], _d[0], _d[1], _d[2] }

// Factory method for Samplers
using samplerCreator = std::function<std::unique_ptr<SamplerBase>(std::string, std::vector<double>)>;
inline GenericFactory< SamplerBase, SamplerTypes, samplerCreator> samplerFactory;

// Base class for all samplers
class SamplerBase{

    protected:
    // Number of samples taken
    int samples = 1;
    // How often should we sample
    int interval;
    // Output filename
    std::string filename;

    public:

    SamplerBase(int interval) : interval(interval){}
    virtual ~SamplerBase() = default;

    virtual void sample(State& state) = 0;
    virtual void save() = 0;
    virtual void close() = 0;
    // Needed by particlesimulator::run
    int getInterval();
};

class Density : public SamplerBase{

    private:
    double binWidth, dh, xb, yb;
    std::vector<unsigned long long int> pDens;
    std::vector<unsigned long long int> nDens;
    int d, bins;
    std::string dim;

    public:

    Density(int d, double dl, double binWidth, double xb, double yb, 
            int interval, std::string filename);
    ~Density() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};

class Energy : public SamplerBase{

    std::vector<double> energies;

    public:
    Energy(int interval, std::string filename);
    ~Energy() override = default;

    void sample(State &state) override;
    void save() override;
    void close() override;
};


class WidomHS : public SamplerBase{
    double cp = 0.0;

    public:

    WidomHS(int interval, std::string filename) : SamplerBase(interval){
        this->filename = "cp_HS_" + filename + ".txt";
    }
    ~WidomHS() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};


class QDist : public SamplerBase{
    private:

    double binWidth;
    std::vector<int> pqDist;
    std::vector<int> nqDist;

    public:

    QDist(double dl, double binWidth, int interval, std::string filename) : SamplerBase(interval){
        this->binWidth = binWidth;
        this->pqDist.resize((int) dl / binWidth, 0);
        this->nqDist.resize((int) dl / binWidth, 0);
        this->filename = filename;
    }
    ~QDist() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};


class XDR : public SamplerBase{
    private:
    XDRFILE *xdf = nullptr;

    public:
    XDR(int interval, std::string filename) : SamplerBase(interval){
        filename = filename + ".xtc";
        xdf = xdrfile_open(filename.c_str(), "w");
    }
    ~XDR() override = default;

    void save() override;
    void close() override;
    void sample(State& state) override;
};


class NumIons : public SamplerBase{
    private:

    std::vector<int> pNum;
    std::vector<int> nNum;

    public:

    NumIons(int interval, std::string filename) : SamplerBase(interval){
        this->filename = filename;
    }
    ~NumIons() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};


class Pressure : public SamplerBase{
    private:

    Eigen::Matrix3d pressureT = Eigen::Matrix3d::Zero();
    double idP = 0.0, V, lZ;
    int counter = 0;
    std::vector<double> pressures;
    public:

    Pressure(int interval, double V, double lZ, std::string filename) : SamplerBase(interval){
        this->filename = filename;
        this->V = V;
        this->lZ = lZ;
        this->filename = "ptZ_" + filename;
        pressures.resize(10000 / interval + 1, 0);
    }
    ~Pressure() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};


class PressureV : public SamplerBase{
    private:

    double av = 0.0, dV, oldV, newV, dL, oldd, old_d;
    std::vector<double> pressures;
    int counter = 0;

    public:

    PressureV(int interval, double dL, double x, double y, double z, std::string filename) : SamplerBase(interval){
        this->filename = filename;
        this->oldV = x*y*z;
        this->newV = x*y*(z+dL);
        this->dV = this->newV - this->oldV;
        //this->oldL = z;
        this->dL = dL;
        this->filename = "pV_" + filename;
        pressures.resize(10000 / interval + 1, 0);
        printf("\t ds: %lf\n", this->dL);
    }
    ~PressureV() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};



class ForcePressure : public SamplerBase{
    private:

    double leftToRight;
    double idP;
    double V;
    double rightToLeft;
    double lZ;
    int counter;
    std::vector<double> leftForce;
    std::vector<double> rightForce;

    public:

    ForcePressure(int interval, double V, double lZ, std::string filename) : SamplerBase(interval){
        this->filename = "fp_" + filename;
        this->V = V;
        this->leftForce.resize(10000 / interval + 1, 0);
        this->rightForce.resize(10000 / interval + 1, 0);
        this->counter = 0;
        this->lZ = lZ;
        this->idP = 0.0;
    }
    ~ForcePressure() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};




class Force : public SamplerBase{
    private:

    Eigen::Vector3d force;

    public:

    Force(int interval, std::string filename) : SamplerBase(interval){
        this->filename = "force_" + filename;
    }
    ~Force() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};





class CliffPressure : public SamplerBase{
    private:
    double dL;
    double rP = 0.0;
    double lP = 0.0;
    double area;
    std::vector<double> pressures;
    int counter = 0;

    public:

    CliffPressure(int interval, double dL, double xL, double yL, std::string filename) : SamplerBase(interval){
        this->filename = "CliffP_" + filename;
        this->dL = dL;
        this->area = xL * yL;
        this->pressures.resize(10000 / interval + 1, 0);
    }
    ~CliffPressure() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};




class ModifiedWidom: public SamplerBase{
    private:
    std::vector<double> nomP{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> denomP{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double facP = 0.0;

    std::vector<double> nomN{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> denomN{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double facN = 0.0;

    int pSamples = 0;
    int nSamples = 0;
    int elecStart = 3;

    public:

    ModifiedWidom(int interval, std::string filename) : SamplerBase(interval){
        this->filename = filename;

        printf("Modified widom sampler");
    }
    ~ModifiedWidom() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};




class ModifiedWidomCoulomb: public SamplerBase{
    private:
    std::vector<double> nomP{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> denomP{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double facP = 0.0;

    std::vector<double> nomN{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> denomN{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double facN = 0.0;

    int pSamples = 0;
    int nSamples = 0;
    int elecStart = 3;

    std::shared_ptr<EnergyBase> e;
    public:

    ModifiedWidomCoulomb(int interval, std::string filename) : SamplerBase(interval){
        this->filename = filename;
        this->e = std::make_shared< PairEnergy<Potentials::Coulomb> >();
        
        e->set_cutoff(120.0);

        printf("Modified widom sampler");
    }
    ~ModifiedWidomCoulomb() override = default;

    void sample(State& state) override;
    void save() override;
    void close() override;
};

} // end of namespace Samplers

} // end of namespace Simulator