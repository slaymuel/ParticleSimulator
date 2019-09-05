#include "particles.h"

class Sampler{
    public:
    int samples;
    virtual void sample(Particles& particles) = 0;
    virtual void save(std::string filename) = 0;
};



class Density : public Sampler{
    private:

    double binWidth, dh;
    std::vector<int> pDens;
    std::vector<int> nDens;
    int d, bins, pSample = 0, nSample = 0;
    public:

    Density(int d, double dl, double binWidth){
        this->binWidth = binWidth;
        this->bins = dl / binWidth + 1;
        this->pDens.resize(this->bins);
        this->nDens.resize(this->bins);
        this->d = d;
        this->dh = dl / 2.0;
    }

    void sample(Particles& particles){
        for(int i = 0; i < particles.tot; i++){
            //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
            if(particles.particles[i]->q > 0){
                pDens[ (int) ((particles.particles[i]->pos[d] + this->dh) / this->binWidth) ]++;
                this->pSample++;
            }
            else{
                nDens[ (int) ((particles.particles[i]->pos[d] + this->dh) / this->binWidth) ]++;
                this->nSample++;
            }
        }

    }

    void save(std::string filename){
        std::ofstream f ("p_" + filename);
        if (f.is_open())
        {
            for(int i = 0; i < this->pDens.size(); i++){
                f << std::fixed << std::setprecision(3) << i * this->binWidth << " " <<  (double) this->pDens[i] / this->pSample << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
        
        std::ofstream fi ("n_" + filename);
        if (fi.is_open())
        {
            for(int i = 0; i < this->nDens.size(); i++){
                fi << std::fixed << std::setprecision(3) << i * this->binWidth << " " <<  (double) this->nDens[i] / this->nSample << "\n";
            }
            fi.close();
        }
        else std::cout << "Unable to open file";
    }
};

