#include "state.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"

class Sampler{

    public:
    int samples = 0;
    int interval;
    std::string filename;

    Sampler(int interval) : interval(interval){}

    virtual void sample(State& state) = 0;
    virtual void save() = 0;
    virtual void close() = 0;
};


namespace Samplers{


class Density : public Sampler{
    private:

    double binWidth, dh, xb, yb;
    std::vector<unsigned long long int> pDens;
    std::vector<unsigned long long int> nDens;
    int d, bins;

    public:

    Density(int d, double dl, double binWidth, double xb, double yb, int interval, std::string filename) : Sampler(interval){
        this->binWidth = binWidth;
        this->bins = dl / binWidth;
        this->pDens.resize(this->bins, 0);
        this->nDens.resize(this->bins, 0);
        this->d = d;    //Which dimension to sample
        this->dh = dl / 2.0;
        this->xb = xb;
        this->yb = yb;
        this->filename = filename;
    }

    void sample(State& state){
        for(unsigned int i = 0; i < state.particles.tot; i++){
            //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
            if(state.particles.particles[i]->q > 0){
                pDens.at( (unsigned int) ( (state.particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
            }

            else{
                nDens.at( (unsigned int) ( (state.particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
            }
        }
        this->samples++;
    }

    void save(){
        std::ofstream f ("p_" + this->filename + ".txt");
        if (f.is_open())
        {
            for(unsigned int i = 0; i < this->pDens.size(); i++){
                f << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 -  this->dh<< " " <<  
                     (double) this->pDens[i] / (this->xb * this->yb * this->binWidth * this->samples) << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
        
        std::ofstream fi ("n_" + this->filename + ".txt");
        if (fi.is_open())
        {
            for(unsigned int i = 0; i < this->nDens.size(); i++){
                fi << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 -  this->dh << " " <<  
                      (double) this->nDens[i] / (this->xb * this->yb * this->binWidth * this->samples) << "\n";
            }
            fi.close();
        }
        else std::cout << "Unable to open file";
    }

    void close(){};
};


class Energy : public Sampler{

    std::vector<double> energies;

    public:
    Energy(int interval, std::string filename) : Sampler(interval){
        energies.reserve(50000);
        this->filename = "energies_" + filename + ".txt";
    }

    void sample(State &state){
        energies.push_back(state.cummulativeEnergy);
    }

    void save(){
        std::ofstream f (this->filename);
        if (f.is_open()){
            for(auto e : energies){
                f << std::fixed << std::setprecision(15) << e << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
    }

    void close(){};
};


class WidomHS : public Sampler{
    double cp = 0.0;

    public:

    WidomHS(int interval, std::string filename) : Sampler(interval){
        this->filename = "cp_HS_" + filename + ".txt";
    }

    void sample(State& state){
        Eigen::Vector3d com = state.geo->random_pos(2.5);
        Eigen::Vector3d qDisp;
        qDisp << 0.0, 0.0, 0.0;
        com[2] = (Random::get_random() * 0.2 - 0.1) * state.geo->_dh[2];
        //std::cout << com[0] << " " << com[1] << " " << com[2] << std::endl;
        state.particles.add(com, com, qDisp, 2.5, state.particles.pModel.rf, state.particles.pModel.q, state.particles.pModel.b, 0.0, 0.0, "WIDOM_PARTICLE");

        if(!state.overlap(state.particles.tot - 1)){
            this->cp += 1.0;
        }
        state.particles.remove(state.particles.tot - 1);
        this->samples++;
        //printf("HS-CP: %lf\n", -std::log(this->cp / this->samples));
        //printf("samples: %d\n", this->samples);
        //printf("cp: %lf\n\n", this->cp);
    }

    void save(){
        std::ofstream f (this->filename);
        if (f.is_open()){
            f << std::fixed << std::setprecision(10) << "Hard-sphere chemical potential: " << -std::log(this->cp / this->samples)  << "\n";
            f.close();
        }
        else std::cout << "Unable to open file";
    }

    void close(){};
};


class QDist : public Sampler{
    private:

    double binWidth;
    std::vector<int> pqDist;
    std::vector<int> nqDist;

    public:

    QDist(double dl, double binWidth, int interval, std::string filename) : Sampler(interval){
        this->binWidth = binWidth;
        this->pqDist.resize((int) dl / binWidth, 0);
        this->nqDist.resize((int) dl / binWidth, 0);
        this->filename = filename;
    }

    void sample(State& state){
        for(unsigned int i = 0; i < state.particles.tot; i++){
            if(state.particles.particles[i]->q > 0.0){

            //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
                //printf("%lf\n", state.geo->distance(state.particles.particles[i]->pos, state.particles.particles[i]->com));
                pqDist.at( (int) ( (state.geo->distance(state.particles.particles[i]->pos, state.particles.particles[i]->com)) /
                                this->binWidth ) )++;
            }
            else{
                nqDist.at( (int) ( (state.geo->distance(state.particles.particles[i]->pos, state.particles.particles[i]->com)) /
                                this->binWidth ) )++; 
            }
        }
        this->samples++;
    }

    void save(){
        std::ofstream pf ("pqDist_" + this->filename + ".txt");
        if (pf.is_open())
        {
            for(unsigned int i = 0; i < this->pqDist.size(); i++){
                pf << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                     (double) this->pqDist[i] / this->samples << "\n";
            }
            pf.close();
        }
        else std::cout << "Unable to open file";

        std::ofstream nf ("nqDist_" + this->filename + ".txt");
        if (nf.is_open())
        {
            for(unsigned int i = 0; i < this->nqDist.size(); i++){
                nf << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                     (double) this->nqDist[i] / this->samples << "\n";
            }
            nf.close();
        }
        else std::cout << "Unable to open file";
    }

    void close(){};
};


class XDR : public Sampler{
    private:
    XDRFILE *xdf = nullptr;

    public:
    XDR(int interval, std::string filename) : Sampler(interval){
        filename = filename + ".xtc";
        xdf = xdrfile_open(filename.c_str(), "w");
    }

    void save(){};

    void close(){
        xdrfile_close(xdf);
    }

    void sample(State& state){
        matrix box;
        box[0][0] = state.geo->_d[0];
        box[0][1] = 0.0;
        box[0][2] = 0.0;
        box[1][0] = state.geo->_d[1];
        box[1][1] = 0.0;
        box[1][2] = 0.0;
        box[2][0] = state.geo->_d[2];
        box[2][1] = 0.0;
        box[2][2] = 0.0;

        
        if (xdf != nullptr) {
            rvec *ps = new rvec[state.particles.tot];
            //size_t N = 0;

            for (unsigned int i = 0; i < state.particles.tot; i++) {
                ps[i][0] = state.particles[i]->pos[0] * 0.1 + state.geo->_d[0] * 0.5;
                ps[i][1] = state.particles[i]->pos[1] * 0.1 + state.geo->_d[1] * 0.5;
                ps[i][2] = state.particles[i]->pos[2] * 0.1 + state.geo->_d[2] * 0.5; 
            }

            write_xtc(xdf, state.particles.tot, state.step, state.step, box, ps, 1000);

            delete[] ps;
        }
        else{
            printf("Could not open xtc file!\n");
            exit(0);
        }
    }
};


class NumIons : public Sampler{
    private:

    std::vector<unsigned long long int> pNum;
    std::vector<unsigned long long int> nNum;

    public:

    NumIons(int interval, std::string filename) : Sampler(interval){
        this->filename = filename;
    }

    void sample(State& state){
        pNum.push_back(state.particles.cTot);
        nNum.push_back(state.particles.aTot);
    }

    void save(){
        std::ofstream f ("pNum_" + this->filename + ".txt");
        if (f.is_open())
        {
            for(unsigned int i = 0; i < this->pNum.size(); i++){
                f << this->pNum[i] << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
        
        std::ofstream fi ("nNum_" + this->filename + ".txt");
        if (fi.is_open())
        {
            for(unsigned int i = 0; i < this->nNum.size(); i++){
                fi << this->nNum[i] << "\n";
            }
            fi.close();
        }
        else std::cout << "Unable to open file";
    }

    void close(){};
};



}

/*
class Potential : public Sampler{
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