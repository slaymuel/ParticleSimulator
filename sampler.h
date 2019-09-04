#include "particles.h"

class Sampler{
    public:
    int samples;
    virtual void sample(Particles& particles) = 0;
    virtual void save(std::string filename) = 0;
};



class Density : public Sampler{
    private:

    double binWidth;
    int bins;
    std::vector<int> density;
    int d;
    double dh;
    public:

    Density(int d, double dl, double binWidth){
        this->binWidth = binWidth;
        this->bins = dl / binWidth + 1;
        this->density.resize(this->bins);
        this->d = d;
        this->dh = dl / 2.0;
    }

    void sample(Particles& particles){
        for(int i = 0; i < particles.tot; i++){
            //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
            density[ (int) ((particles.particles[i]->pos[d] + this->dh) / this->binWidth) ]++;
        }

        this->samples++;
    }

    void save(std::string filename){
        std::ofstream f (filename);
        if (f.is_open())
        {
            for(int i = 0; i < this->density.size(); i++){
                f << std::fixed << std::setprecision(3) << i * this->binWidth << " " <<  (double) this->density[i] / this->samples << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
    }
};

