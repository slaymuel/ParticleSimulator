#include "state.h"

class Sampler{

    public:
    int samples = 0;
    int interval;

    Sampler(int interval) : interval(interval){}

    virtual void sample(State& state) = 0;
    virtual void save(std::string filename) = 0;
};


namespace Samplers{


class Density : public Sampler{
    private:

    double binWidth, dh, xb, yb;
    std::vector<int> pDens;
    std::vector<int> nDens;
    int d, bins;

    public:

    Density(int d, double dl, double binWidth, double xb, double yb, int interval) : Sampler(interval){
        this->binWidth = binWidth;
        this->bins = dl / binWidth;
        this->pDens.resize(this->bins, 0);
        this->nDens.resize(this->bins, 0);
        this->d = d;    //Which dimension to sample
        this->dh = dl / 2.0;
        this->xb = xb;
        this->yb = yb;
    }

    void sample(State& state){
        for(unsigned int i = 0; i < state.particles.tot; i++){
            //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
            if(state.particles.particles[i]->q > 0){
                pDens.at( (int) ( (state.particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
            }

            else{
                nDens.at( (int) ( (state.particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
            }
        }
        this->samples++;
    }

    void save(std::string filename){
        std::ofstream f ("p_" + filename + ".txt");
        if (f.is_open())
        {
            for(unsigned int i = 0; i < this->pDens.size(); i++){
                f << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                     (double) this->pDens[i] / (this->xb * this->yb * this->binWidth * this->samples) << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
        
        std::ofstream fi ("n_" + filename + ".txt");
        if (fi.is_open())
        {
            for(unsigned int i = 0; i < this->nDens.size(); i++){
                fi << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                      (double) this->nDens[i] / (this->xb * this->yb * this->binWidth * this->samples) << "\n";
            }
            fi.close();
        }
        else std::cout << "Unable to open file";
    }
};


class Energy : public Sampler{

    std::vector<double> energies;

    public:
    Energy(int interval) : Sampler(interval){
        energies.reserve(50000);
    }

    void sample(State &state){
        energies.push_back(state.cummulativeEnergy);
    }

    void save(std::string filename){
        std::ofstream f ("energies_" + filename + ".txt");
        if (f.is_open()){
            for(auto e : energies){
                f << std::fixed << std::setprecision(15) << e << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
    }
};


class WidomHS : public Sampler{
    double cp = 0.0;

    public:

    WidomHS(int interval) : Sampler(interval){}

    void sample(State& state){
        Eigen::Vector3d com = state.geo->random_pos(2.5);
        com[2] = (Random::get_random() * 0.2 - 0.1) * state.geo->dh[2];
        //std::cout << com[0] << " " << com[1] << " " << com[2] << std::endl;
        state.particles.add(com, com, 2.5, state.particles.pModel.rf, state.particles.pModel.q, state.particles.pModel.b, "WIDOM_PARTICLE");
        if(!state.overlap(state.particles.tot - 1)){
            this->cp += 1.0;
        }
        state.particles.remove(state.particles.tot - 1);
        this->samples++;
        printf("HS-CP: %lf\n", -std::log(this->cp / this->samples));
        printf("samples: %d\n", this->samples);
        printf("cp: %lf\n\n", this->cp);
    }

    void save(std::string filename){
        std::ofstream f ("cp_HS_" + filename + ".txt");
        if (f.is_open()){
            f << std::fixed << std::setprecision(10) << "Hard-sphere chemical potential: " << -std::log(this->cp / this->samples)  << "\n";
            f.close();
        }
        else std::cout << "Unable to open file";
    }
};


class QDist : public Sampler{
    private:

    double binWidth;
    std::vector<int> pqDist;
    std::vector<int> nqDist;

    public:

    QDist(double dl, double binWidth, int interval) : Sampler(interval){
        this->binWidth = binWidth;
        this->pqDist.resize((int) dl / binWidth, 0);
        this->nqDist.resize((int) dl / binWidth, 0);
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

    void save(std::string filename){
        std::ofstream pf ("pqDist_" + filename + ".txt");
        if (pf.is_open())
        {
            for(unsigned int i = 0; i < this->pqDist.size(); i++){
                pf << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                     (double) this->pqDist[i] / this->samples << "\n";
            }
            pf.close();
        }
        else std::cout << "Unable to open file";

        std::ofstream nf ("nqDist_" + filename + ".txt");
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