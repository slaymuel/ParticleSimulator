class Sampler{
    int samples;
};



class Density : public Sampler{
    double binWidth;
    int bins;
    std::vector<double> density;

    Density(int d, double dl, double binWidth){
        this->binWidth = binWidth;
        this->bins = dl / binWidth + 1;
        this->density.resize(this->bins);
    }

    void sample(){
        for(){

        }

        this->samples++;
    }

    void save(std::string filename){

    }
};

