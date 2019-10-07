#include "particles.h"

class Sampler{
    public:
    int samples = 0;
    virtual void sample(Particles& particles) = 0;
    virtual void save(std::string filename) = 0;
};



class Density : public Sampler{
    private:

    double binWidth, dh, xb, yb;
    std::vector<int> pDens;
    std::vector<int> nDens;
    int d, bins;
    public:

    Density(int d, double dl, double binWidth, double xb, double yb){
        this->binWidth = binWidth;
        this->bins = dl / binWidth + 1;
        this->pDens.resize(this->bins);
        this->nDens.resize(this->bins);
        this->d = d;    //Which dimension to sample
        this->dh = dl / 2.0;
        this->xb = xb;
        this->yb = yb;
    }

    void sample(Particles& particles){
        for(unsigned int i = 0; i < particles.tot; i++){
            //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
            if(particles.particles[i]->q > 0){
                pDens.at( (int) ( (particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
            }

            else{
                nDens.at( (int) ( (particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
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

