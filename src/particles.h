#pragma once

#include <tuple>
#include <map>
#include "particle.h"
#include "logger.h"

namespace Simulator{

// Particles object holds all particles
class Particles{

    private:
    using ParticleVec = std::vector< std::shared_ptr<Particle> >;

    // If the models are created
    bool setPModel = false;
    bool setNModel = false;
    ParticleVec particles;

    public:

    // Not used yet
    Eigen::MatrixXd positions;
    // Model particles
    Particle pModel;
    Particle nModel;
    // Optimize by setting fixed capacity of particles vector
    // Holds moved particles after each move
    std::vector<int> movedParticles;
    //The number of particles
    unsigned int cTot = 0, aTot = 0, tot = 0;

    Particles();

    // Iterator behaviour for range based loops etc
    ParticleVec::iterator begin();
    ParticleVec::iterator end();
    ParticleVec::const_iterator begin() const;
    ParticleVec::const_iterator end() const;

    bool is_valid();
    // Get a subset of particle positions using positions matrix
    Eigen::MatrixXd get_subset(int sr, int fr) const;
    // Get position of a particles using positions matrix
    Eigen::MatrixXd get_particle_pos(int index) const;
    // Indexing operators
    std::shared_ptr<Particle> operator[](std::size_t index);
    const std::shared_ptr<Particle> operator[](std::size_t index) const;
    // Get the length of the particles vector
    unsigned int size() const;
    // Get a subset of particles
    std::vector< std::shared_ptr<Particle> > get_subset(std::vector<unsigned int> &ps) const;
    std::vector< std::shared_ptr<Particle> > get_subset(unsigned int i) const;
    // Get a random particle
    std::shared_ptr<Particle> random() const;
    // Translate a group of particles
    void translate(std::vector<unsigned int> &ps, std::vector<double> &disp);
    void translate(std::vector<unsigned int> &ps, Eigen::Vector3d &disp);
    // Remove a random particle
    std::tuple<unsigned int, double> remove_random();
    // Remove a random particle of specie "type"
    std::tuple<unsigned int, double> remove_random(double type);
    // Remove particle at index "index"
    void remove(std::size_t index);
    // Create pNum cations and nNum anions with parameters params
    void create(int pNum, int nNum, std::map<std::string, double> params);
    // Add a random particle specie to a random position in the box
    // Returns the index and charge of the added particle
    std::tuple<unsigned int, double> add_random(std::vector<double> box, double type = 0.0);
    // Add specific particle
    void add(const std::shared_ptr<Particle>& p, int index = -1);
    // Add particle with specified parameters
    template <typename T, typename G>
    void add(T com, T pos, G qDisp, double r, double rf, double q, double b, double b_min, double b_max, std::string name){
        this->particles.push_back(std::make_shared<Particle>());

        this->particles[this->tot]->com << com[0], com[1], com[2];
        this->particles[this->tot]->pos << pos[0], pos[1], pos[2];
        this->particles[this->tot]->qDisp << qDisp[0], qDisp[1], qDisp[2];

        if(std::abs(this->particles[this->tot]->qDisp.norm() - b) > 1e-10){
            printf("Error reading particles, b is not equal to |qDisp| for particle %u\n", this->tot);
            printf("b = %.10lf\n", b);
            printf("|qDisp| = %.10lf\n", this->particles[this->tot]->qDisp.norm());
            exit(1);   
        }

        this->particles[this->tot]->index = this->tot;
        this->particles[this->tot]->r = r;
        this->particles[this->tot]->rf = rf;
        this->particles[this->tot]->q = q;
        this->particles[this->tot]->b = b;
        this->particles[this->tot]->b_min = b_min;
        this->particles[this->tot]->b_max = b_max;
        this->particles[this->tot]->name = name;
            
        if(q > 0)
            this->cTot++;
        else 
            this->aTot++;

        this->tot++;
        assert(this->tot <= this->particles.size() && "tot is larger than particle vector size\n");

        if(!this->setPModel){
            if(q > 0){
                this->pModel.q = q;
                this->pModel.b_min = b_min;
                this->pModel.b_max = b_max;
                this->pModel.r = r;
                this->pModel.rf = rf;
                this->pModel.name = name;
                this->setPModel = true;
            }
        }

        if(!this->setNModel){
            if(q < 0){
                this->nModel.q = q;
                this->nModel.b_min = b_min;
                this->nModel.b_max = b_max;
                this->nModel.r = r;
                this->nModel.rf = rf;
                this->nModel.name = name;
                this->setNModel = true;
            }
        }
    }


    template <typename T>
    void add(T com, double r, double rf, double q, double b_min, double b_max, std::string name){
        Eigen::Vector3d qDisp = Random::get_random_vector(b_max);
        Eigen::Vector3d pos = com + qDisp;
        double b = qDisp.norm();
        add(com, pos, qDisp, r, rf, q, b, b_min, b_max, name);
    }
};

}