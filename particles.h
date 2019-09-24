#pragma once

#include <Eigen/Dense>
#include <vector>
#include "particle.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

class Particles{

    public:
    //Eigen::MatrixXd positions;
    Particle pModel;
    Particle nModel;
    std::vector< std::shared_ptr<Particle> > particles, cations, anions;
    std::vector<int> movedParticles;
    unsigned int pTot = 0, cTot = 0, aTot = 0, iTot = 0, tot = 0;

    //Eigen::MatrixXd get_subset(int sr, int fr){
    //    return this->positions.block(sr, 0, fr, 3);
    //}

    Particles(){}

    std::shared_ptr<Particle> operator[](std::size_t index){
        return particles[index];
    }

    std::vector<std::shared_ptr<Particle>> get_subset(std::vector<unsigned int> &ps){
        std::vector<std::shared_ptr<Particle>> subset;
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


    template <typename T>
    void add(T pos, double r, double rf, double q, double b, std::string name, bool image = false){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];

        //Create particle
        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[this->tot]->pos  << pos[0], pos[1], pos[2];
        //std::cout << this->particles[tot]->pos << " " << std::endl;
        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->name = name;


        if(!image){
            
            if(q > 0){
                this->cTot++;
            }
            else {
                this->aTot++;
            }

            this->pTot++;
        }
        else{
            this->iTot++;
        }

        this->tot++;
        assert(tot <= this->particles.size() && "tot is larger than particle vector size\n");
        assert(tot == this->pTot + this->iTot && "pTot + iTot != tot!\n");
    }


    void add(const std::shared_ptr<Particle>& p){
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];
        //Create particle
        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }
        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[tot]->pos  << p->pos[0], p->pos[1], p->pos[2];
        //std::cout << this->particles[tot]->pos << " " << std::endl;
        this->particles[tot]->index = this->tot;
        this->particles[tot]->r = p->r;
        this->particles[tot]->rf = p->rf;
        this->particles[tot]->q = p->q;
        this->particles[tot]->b = p->b;
        this->particles[tot]->name = p->name;

        if(p->q > 0){
            this->cTot++;
        }
        else {
            this->aTot++;
        }

        this->pTot++;
        this->tot++;
    }


    void add(const std::shared_ptr<Particle>& p, int index){
        
        //Resize positions
        //this->positions.conservativeResize(this->positions.rows() + 1, 3);
        //this->positions.row(this->positions.rows() - 1) << pos[0], pos[1], pos[2];
        //Create particle

        if(this->tot == this->particles.size()){
            //printf("Allocating\n");
            this->particles.push_back(std::make_shared<Particle>());
        }
        //std::copy(this->particles.begin() + index, this->particles.begin() + this->tot, this->particles.begin() + index + 1);
        for(int i = this->tot; i > index; i--){
            //printf("Moving particle %i to %i\n", i - 1, i);
            *(this->particles[i]) = *(this->particles[i - 1]);
            this->particles[i]->index = i;
        }


        this->tot++;

        //this->particles.back()->pos = this->positions.row(this->positions.rows() - 1);
        this->particles[index]->pos << p->pos[0], p->pos[1], p->pos[2];
        //std::cout << this->particles[tot]->pos << " " << std::endl;
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

        this->pTot++;
    }


    void remove(std::size_t index){
        //this->particles.erase(this->particles.begin() + i);
        //printf("Copying in remove\n");
        //std::copy(this->particles.begin() + i + 1, this->particles.begin() + this->tot, this->particles.begin() + i);
        //std::copy(this->particles.begin() + i, this->particles.begin() + this->tot, this->particles.begin() + i - 1);

        //this->particles.erase(this->particles.begin() + i);
        if(this->particles[index]->q > 0){
            this->cTot--;
        }
        else{
            this->aTot--; 
        }


        for(unsigned int i = index; i < this->tot - 1; i++){
            //printf("Remove: Moving particle %i to %i\n", i + 1, i);
            *(this->particles[i]) = *(this->particles[i + 1]);
            this->particles[i]->index = i;
        }

        this->pTot--;
        this->tot--;


    }



    void load(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<std::string> names){
        //assert correct sizes
        bool setNModel = false, setPModel = false;
        for(unsigned int i = 0; i < pos.size(); i++){

            this->add(pos[i], r[i], rf[i], charges[i], b[i], names[i]);
            
            if(!setPModel){
                if(charges[i] > 0){
                    pModel.q = charges[i];
                    pModel.r = r[i];
                    pModel.rf = rf[i];
                    setPModel = true;
                }
            }
            if(!setNModel){
                if(charges[i] < 0){
                    nModel.q = charges[i];
                    nModel.r = r[i];
                    nModel.rf = rf[i];
                    setNModel = true;
                }
            }

        }
        if(!setPModel || !setNModel){
            printf("Cation or anion model not set!\n");
            exit(1);
        }
        printf("Loaded %u particles\n", this->tot);
    }



    void create(int pNum, int nNum, double p, double n, double rfp = 2.5, double rfn = 2.5){
        pModel.q = p;
        pModel.r = 2.5;
        pModel.rf = rfp;

        nModel.q = n;
        nModel.r = 2.5;
        nModel.rf = rfn;

        std::vector<double> v;
        for(int i = 0; i < pNum + nNum; i++){
            v = Random::get_vector();
            (i < pNum) ? this->add(v, 2.5, rfp, p, 0.0, "Na") : this->add(v, 2.5, rfn, n, 0.0, "Cl") ;
        }
        printf("\nCreated %i cations and %i anions\n", pNum, nNum);
    }


    //Write xyz file
    void to_xyz(std::string fileName){
        std::ofstream f (fileName + ".xyz");
        if (f.is_open())
        {
            f << this->tot << "\n\n";
            for(unsigned int i = 0; i < this->tot; i++){

                f << std::fixed << std::setprecision(3) << this->particles[i]->name << " " <<  this->particles[i]->pos[0] << " " << this->particles[i]->pos[1] << " " << this->particles[i]->pos[2] << "\n";
            }
            f << "10 10 10" << "\n";
            f.close();
        }
        else std::cout << "Unable to open file";
    }

    //write checkpoint file
    void to_cpt(std::string fileName){
        std::ofstream f (fileName + ".cp");
        if (f.is_open())
        {
            for(unsigned int i = 0; i < this->tot; i++){

                f << std::fixed << std::setprecision(15) << " " <<  this->particles[i]->pos[0] << " " << this->particles[i]->pos[1] << " " << this->particles[i]->pos[2] << " " << 
                                                                    this->particles[i]->pos[0] << " " << this->particles[i]->pos[1] << " " << this->particles[i]->pos[2] << " " << 
                                                                    this->particles[i]->q << " " << this->particles[i]->r << " " << this->particles[i]->rf << " " << 
                                                                    this->particles[i]->b << " " << this->particles[i]->name << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
    }
};