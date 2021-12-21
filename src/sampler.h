#pragma once

#include "state.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"

namespace Samplers{

class SamplerBase{

    public:
    int samples = 1;
    int interval;
    std::string filename;

    SamplerBase(int interval) : interval(interval){}

    virtual void sample(State& state) = 0;
    virtual void save() = 0;
    virtual void close() = 0;
};


class Density : public SamplerBase{
    private:

    double binWidth, dh, xb, yb;
    std::vector<unsigned long long int> pDens;
    std::vector<unsigned long long int> nDens;
    int d, bins;
    std::string dim;

    public:

    Density(int d, double dl, double binWidth, double xb, double yb, int interval, std::string filename) : SamplerBase(interval){
        this->binWidth = binWidth;
        this->bins = dl / binWidth;
        this->pDens.resize(this->bins, 0);
        this->nDens.resize(this->bins, 0);
        this->d = d;    //Which dimension to sample
        this->dh = dl / 2.0;
        this->xb = xb;
        this->yb = yb;
        this->filename = filename;
        if(d == 0) this->dim = "x";
        if(d == 1) this->dim = "y";
        if(d == 2) this->dim = "z";
    }

    void sample(State& state);
    void save();
    void close();
};


class Energy : public SamplerBase{

    std::vector<double> energies;

    public:
    Energy(int interval, std::string filename) : SamplerBase(interval){
        energies.reserve(50000);
        this->filename = "energies_" + filename + ".txt";
    }

    void sample(State &state);
    void save();
    void close();
};


class WidomHS : public SamplerBase{
    double cp = 0.0;

    public:

    WidomHS(int interval, std::string filename) : SamplerBase(interval){
        this->filename = "cp_HS_" + filename + ".txt";
    }

    void sample(State& state);
    void save();
    void close();
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

    void sample(State& state);
    void save();
    void close();
};


class XDR : public SamplerBase{
    private:
    XDRFILE *xdf = nullptr;

    public:
    XDR(int interval, std::string filename) : SamplerBase(interval){
        filename = filename + ".xtc";
        xdf = xdrfile_open(filename.c_str(), "w");
    }

    void save();
    void close();
    void sample(State& state);
};


class NumIons : public SamplerBase{
    private:

    std::vector<int> pNum;
    std::vector<int> nNum;

    public:

    NumIons(int interval, std::string filename) : SamplerBase(interval){
        this->filename = filename;
    }

    void sample(State& state);
    void save();
    void close();
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

    void sample(State& state);
    void save();
    void close();
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

    void sample(State& state);
    void save();
    void close();
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

    void sample(State& state);
    void save();
    void close();
};







class Force : public SamplerBase{
    private:

    Eigen::Vector3d force;

    public:

    Force(int interval, std::string filename) : SamplerBase(interval){
        this->filename = "force_" + filename;
    }

    void sample(State& state);
    void save();
    void close();
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

    void sample(State& state);
    void save();
    void close();
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

    void sample(State& state);
    void save();
    void close();
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
        this->e = std::make_shared< PairEnergy<Coulomb> >();
        
        e->set_cutoff(120.0);

        printf("Modified widom sampler");
    }

    void sample(State& state);
    void save();
    void close();
};











}

/*
class Potential : public SamplerBase{
    private:

    double binWidth, dl;
    std::vector< std::vector<double> > potential;
    int bins;

    public:

    Potential(double dl, double binWidth){
        this->binWidth = binWidth;
        this->bins = dl / binWidth;
        this->potential.resize(this->bins, std::vector<double>(2));
        this->dl = dl;

        for(auto p : this->potential){
            p[0] = 0.0;
            p[1] = 0.0;
        }
    }

    void sample(Particles& particles){
        particles.add(T com, T pos, double r, 2.5, 1.0, 0.0, std::string name, bool image = false);
        this->potential[(int)particles.particles.back()->pos[2] + 0.5 * this->dl][0] += get_energy();
        this->potential[(int)particles.particles.back()->pos[2] + 0.5 * this->dl][1] += 1;
        particles.remove(particles.particles.back()->index);
    }

    void save(std::string filename){
        std::ofstream f ("potential_" + filename + ".txt");
        if (f.is_open())
        {
            for(unsigned int i = 0; i < this->potential.size(); i++){
                f << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                     (double) this->potential[i][0] / this->potential[i][1] << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
    }
};
*/