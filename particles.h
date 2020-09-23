#pragma once

#include <Eigen/Dense>
#include <vector>
#include "particle.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <tuple>
//#include "../libxdrfile/include/xdrfile_xtc.h"

class Particles{

    public:
    //Eigen::MatrixXd positions;
    Particle pModel;
    Particle nModel;
    std::vector< std::shared_ptr<Particle> > particles, cations, anions;
    std::vector<int> movedParticles;
    unsigned int cTot = 0, aTot = 0, tot = 0;

    //Eigen::MatrixXd get_subset(int sr, int fr){
    //    return this->positions.block(sr, 0, fr, 3);
    //}

    Particles(){}

    std::shared_ptr<Particle> operator[](std::size_t index){
        return particles[index];
    }

    std::vector< std::shared_ptr<Particle> > get_subset(std::vector<unsigned int> &ps){
        std::vector< std::shared_ptr<Particle> > subset;

        //std::vector<Particle> subset(ps.size(), 0);
        //std::transform(ps.begin(), ps.end(), subset.begin(), [particles](size_t i) {return particles[i];});
        for (auto p : ps){
            if (p < this->tot) subset.push_back(this->particles[p]);
        }
        return subset;
    }

    std::shared_ptr<Particle> random(){
        //https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
        //std::sample

        //return this->particles[(*distribution)(rand_gen)];
        return this->particles[Random::get_random(this->tot)];
    }

    void translate(std::vector<unsigned int> &ps, std::vector<double> &disp){
        for(auto p : ps){
            this->particles[p]->translate(disp);
        }
    }

    void translate(std::vector<unsigned int> &ps, Eigen::Vector3d &disp){
        for(auto p : ps){
            this->particles[p]->translate(disp);
        }
    }

    void set_models(std::vector<double> q, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<std::string> names){
        this->pModel.q =        q[0];
        this->pModel.r =        r[0];
        this->pModel.rf =      rf[0];
        this->pModel.b =        b[0];
        this->pModel.name = names[0]; 

        this->nModel.q =        q[1];
        this->nModel.r =        r[1];
        this->nModel.rf =      rf[1];
        this->nModel.b =        b[1];
        this->nModel.name = names[1]; 
    }

    template <typename T>
    void add(T com, T pos, double r, double rf, double q, double b, std::string name){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[this->tot]->com << com[0], com[1], com[2];
        this->particles[this->tot]->pos << pos[0], pos[1], pos[2];
        this->particles[this->tot]->qDisp = this->particles[this->tot]->pos - this->particles[this->tot]->com;
        this->particles[this->tot]->qDisp = this->particles[this->tot]->qDisp.stableNormalized() * b;
        //this->particles[this->tot]->pos = this->particles[this->tot]->com + this->particles[this->tot]->qDisp;

        //std::cout << this->particles[tot]->pos << " " << std::endl;
        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->name = name;


            
        if(q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->tot++;
        assert(this->tot <= this->particles.size() && "tot is larger than particle vector size\n");
    }

    template <typename T>
    void add(T com, double r, double rf, double q, double b, std::string name){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[this->tot]->com << com[0], com[1], com[2];
        this->particles[this->tot]->qDisp = Random::get_norm_vector();
        this->particles[this->tot]->qDisp = this->particles[this->tot]->qDisp.stableNormalized() * b;
        this->particles[this->tot]->pos = this->particles[this->tot]->com + this->particles[this->tot]->qDisp;
        //std::cout << this->particles[tot]->pos << " " << std::endl;
        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->name = name;

        if(q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->tot++;
        assert(this->tot <= this->particles.size() && "tot is larger than particle vector size\n");
    }


    void add(const std::shared_ptr<Particle>& p){

        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}

        this->particles[this->tot]->com << p->com[0], p->com[1], p->com[2];
        this->particles[this->tot]->pos << p->pos[0], p->pos[1], p->pos[2];
        this->particles[this->tot]->qDisp << p->qDisp;

        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = p->r;
        this->particles[this->tot]->rf = p->rf;
        this->particles[this->tot]->q = p->q;
        this->particles[this->tot]->b = p->b;
        this->particles[this->tot]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->tot++;
    }


    void add(const std::shared_ptr<Particle>& p, int index){

        //if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        //}
        //std::copy(this->particles.begin() + index, this->particles.begin() + this->tot, this->particles.begin() + index + 1);
        for(int i = this->tot; i > index; i--){
            //printf("Moving particle %i to %i\n", i - 1, i);
            *(this->particles[i]) = *(this->particles[i - 1]);
            this->particles[i]->index = i;
        }

        this->tot++;

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[index]->com << p->com[0], p->com[1], p->com[2];
        this->particles[index]->pos << p->pos[0], p->pos[1], p->pos[2];
        this->particles[index]->qDisp << p->qDisp;

        this->particles[index]->index = index;
        this->particles[index]->r = p->r;
        this->particles[index]->rf = p->rf;
        this->particles[index]->q = p->q;
        this->particles[index]->b = p->b;
        this->particles[index]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }      
    }



    std::tuple<unsigned int, double> add_random(std::vector<double> box){
        Eigen::Vector3d com;
        double rand = Random::get_random();
        double q;

        //Add cation
        if(rand < 0.5){
            com = Random::random_pos_box(this->pModel.rf, box);
            this->add(com, this->pModel.r, this->pModel.rf, this->pModel.q, this->pModel.b, "Na");
            q = this->pModel.q;
        }

        //Add anion
        else{
            com = Random::random_pos_box(this->nModel.rf, box);
            this->add(com, this->nModel.r, this->nModel.rf, this->nModel.q, this->nModel.b, "Cl");
            q = this->nModel.q;
        }
        //printf("Adding %i charge %lf com %lf %lf %lf\n", this->tot - 1, q, com[0], com[1], com[2]);
        return {this->tot - 1, q};
    }

    std::tuple<unsigned int, double> remove_random(){
        double q;
        double rand = Random::get_random();
        int rand2 = Random::get_random(this->tot);
        //int rand4 = Random::get_random(this->tot);
        //int rand2 = (int) rand3;
        //printf("rand3 %lf rand2 %i rand4 %i\n", rand3, rand2, rand4);
        //printf("rand2 %i\n", rand2);
        if(rand < 0.5){
            if(this->cTot > 0){
                do{
                    rand2 = Random::get_random(this->tot);
                } while(this->particles[rand2]->q != this->pModel.q);
            }
        }
        else{
            if(this->aTot > 0){
                do{
                    rand2 = Random::get_random(this->tot);
                } while(this->particles[rand2]->q != this->nModel.q);
            }
        }

        q = this->particles[rand2]->q;
        this->remove(rand2); // remove particle
        if(this->cTot <= 0 || this->aTot <= 0){
            printf("Woops, no particles left!\n");
            exit(0);
        }
        return {rand2, q};
        //printf("Removing %i charge %lf\n", rand2, q);
        
    }



    void remove(std::size_t index){
        //this->particles.erase(this->particles.begin() + i);
        //printf("Copying in remove\n");
        //std::copy(this->particles.begin() + i + 1, this->particles.begin() + this->tot, this->particles.begin() + i);
        //std::copy(this->particles.begin() + i, this->particles.begin() + this->tot, this->particles.begin() + i - 1);

        if(this->particles[index]->q > 0){
            this->cTot--;
        }
        else{
            this->aTot--; 
        }
        
        //move last particle in particles to position of particle to be removed
        this->particles.erase(this->particles.begin() + index);

        for(unsigned int i = index; i < this->tot - 1; i++){
            //printf("Remove: Moving particle %i to %i\n", i + 1, i);
            //*(this->particles[i]) = *(this->particles[i + 1]);
            this->particles[i]->index = i;
        }

        this->tot--;
    }



    void load(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<std::string> names){
        //assert correct sizes
        bool setNModel = false, setPModel = false;
        for(unsigned int i = 0; i < pos.size(); i++){

            this->add(com[i], pos[i], r[i], rf[i], charges[i], b[i], names[i]);
            //this->add(com[i], com[i], 2.5, rf[i], charges[i], 0.0, names[i]);
            
            if(!setPModel){
                if(charges[i] > 0){
                    pModel.q = charges[i];
                    pModel.b = b[i];
                    pModel.r = r[i];
                    pModel.rf = rf[i];
                    pModel.name = names[i];
                    setPModel = true;
                }
            }

            if(!setNModel){
                if(charges[i] < 0){
                    nModel.q = charges[i];
                    nModel.b = b[i];
                    nModel.r = r[i];
                    nModel.rf = rf[i];
                    nModel.name = names[i];
                    setNModel = true;
                }
            }

        }
        if(!setPModel || !setNModel){
            printf("Cation or anion model not set!\n");
            exit(1);
        }
        printf("Loaded %u particles, %u cations and %u anions.\n", this->tot, this->cTot, this->aTot);
    }



    void create(int pNum, int nNum, double p, double n, double rfp = 2.5, double rfn = 2.5, double rp = 2.5, double rn = 2.5, double bp = 0.0, double bn = 0.0){

        pModel.q = p;
        pModel.r = rp;
        pModel.rf = rfp;
        pModel.b = bp;
        pModel.name = "Na";

        nModel.q = n;
        nModel.r = rn;
        nModel.rf = rfn;
        nModel.b = bn;
        nModel.name = "Cl";

        Eigen::Vector3d com;
        for(int i = 0; i < pNum + nNum; i++){
            com = Random::get_vector();
            (i < pNum) ? this->add(com, rp, rfp, p, bp, "Na") : this->add(com, rn, rfn, n, bn, "Cl");
        }
        printf("\nCreated %i cations and %i anions\n", pNum, nNum);
    }


    void read_cp(std::string fileName){
        std::vector< std::vector<double> > data;
        std::ifstream file(fileName);
        std::string line;
        std::string name;
        int linenr = 0;
        bool setPModel = false, setNModel = false;

        printf("\nLoading checkpoint file \"%s\"....\n", fileName.c_str());
        while(std::getline(file, line)){
            //printf("Getting line: %i\n", linenr);

            std::vector<double>   lineData;
            std::stringstream  lineStream(line);

            double value;
            //printf("Reading values\n");
            while(lineStream >> value){
                //printf("%lf ", value);
                lineData.push_back(value);
            }
            //printf("\n");
//void add(T com, T pos, double r, double rf, double q, double b, std::string name, bool image = false){
            name = (lineData[8] > 0) ? "Na" : "Cl";
            std::vector<double> com{lineData[0], lineData[1], lineData[2]};
            std::vector<double> pos{lineData[3], lineData[4], lineData[5]};
            /*
            lineData[6] = q
            lineData[7] = r
            lineData[8] = rf
            lineData[9] = b
            lineData[10] = name
            */
            //printf("Adding particle\n\n");
            this->add(com, pos, lineData[7], lineData[8], lineData[6], lineData[9], name);

            if(!setPModel){
                if(lineData[6] > 0){
                    pModel.q = lineData[6];
                    pModel.b = lineData[9];
                    pModel.r = lineData[7];
                    pModel.rf = lineData[8];
                    pModel.name = name;
                    setPModel = true;
                }
            }

            if(!setNModel){
                if(lineData[6] < 0){
                    nModel.q = lineData[6];
                    nModel.b = lineData[9];
                    nModel.r = lineData[7];
                    nModel.rf = lineData[8];
                    nModel.name = name;
                    setNModel = true;
                }
            }


            //data.push_back(lineData);
            linenr++;
        }
        printf("Loaded %i particles\n", linenr);
    }
};