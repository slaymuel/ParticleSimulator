#pragma once

#include "particle.h"
#include "particles.h"
#include "geometry.h"


class EnergyBase{

    protected:

    
    double cutoff;

    public:
    Geometry *geo;

    //is not used by ext - bad design
    //should add one layer EnergyBase -> PairEnergy -> Pair
    //where PairEnergy constains set_geo()
    // EnergyBase -> External -> ExtEnergy
    void set_geo(Geometry* geo){
        this->geo = geo;
    }

    void set_cutoff(double cutoff){
        this->cutoff = cutoff;
        printf("\tEnergy cutoff set to %lf\n", this->cutoff);
    }


    virtual ~EnergyBase(){};
    virtual double all2all(Particles& particles) = 0;
    virtual double i2all(std::shared_ptr<Particle> p, Particles& particles) = 0;
    virtual double operator()(std::vector< unsigned int >&& p, Particles& particles) = 0;
    virtual Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles) = 0;
    virtual Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) = 0;
    virtual Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2) = 0;
    virtual double operator()(std::vector< unsigned int >& p, Particles& particles) = 0;
    virtual void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new) = 0;
    virtual void update(double x, double y, double z) = 0;
    virtual void initialize(Particles& particles) = 0;
};



template <typename E>
class PairEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    double all2all(Particles& particles){
        double e = 0.0;

        //printf("all2all geo: %lf %lf %lf\n", this->geo->dh[0], this->geo->dh[1], this->geo->dh[2]);
        #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = i + 1; j < particles.tot; j++){
                //printf("distance %lf\n", this->geo->distance(particles[i]->pos, particles[j]->pos));
                //printf("Indices: %u, %u\n", i, j);
                e += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
            }  
        }
        //printf("Real energy: %.15lf\n", e);
        return e * constants::lB;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = 0.0;
        
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            e += i2i(p->q, particles[i]->q, this->geo->distance(p->pos, particles[i]->pos));
        }

        return e;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        //printf("i2all geo: %lf %lf %lf\n", this->geo->dh[0], this->geo->dh[1], this->geo->dh[2]);
        double e = 0.0;
        // Need to fix this, not a nice solution (when volume move).........
        if(p.size() == particles.tot){
            e = all2all(particles) / constants::lB;
        }
        else{
            #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
            for(int i = 0; i < p.size(); i++){
                //do instead i2all(s, particles);
                e += i2all(particles.particles[p[i]], particles);
            }

            for(int i = 0; i < p.size(); i++){
                for(int j = i + 1; j < p.size(); j++){
                    e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
                }
            }
        }

        return e * constants::lB;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){
        //printf("i2all geo: %lf %lf %lf\n", this->geo->dh[0], this->geo->dh[1], this->geo->dh[2]);
        double e = 0.0;
        // Need to fix this, not a nice solution.........
        if(p.size() == particles.tot){
            e = all2all(particles) / constants::lB;
        }
        else{
            #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
            for(int i = 0; i < p.size(); i++){
                //do instead i2all(s, particles);
                e += i2all(particles.particles[p[i]], particles);
            }

            for(int i = 0; i < p.size(); i++){
                for(int j = i + 1; j < p.size(); j++){
                    e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
                }
            }
        }

        return e * constants::lB;
    }

    inline double i2i(double& q1, double& q2, double&& dist){
        if(dist <= this->cutoff){
            return energy_func(q1, q2, dist);
        }
        else{
            return 0.0;
        }
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void initialize(Particles& particles){}
    void update(double x, double y, double z){}

    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for(unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            force += energy_func.force(p->q, particles[i]->q, this->geo->displacement(p->pos, particles[i]->pos));
        }
        return force * constants::lB;
    }

    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force * constants::lB;
        return force(p1->pos, p2->pos, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d disp = this->geo->displacement(pos1, pos2);
        if(disp.norm() <= this->cutoff){
            return energy_func.force(q1, q2, disp) * constants::lB;
        }
        else{
            return Eigen::Vector3d::Zero();
        }

    }
};









template <typename E>
class PairEnergyCOM : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:
    PairEnergyCOM(double k, double R){
        energy_func.set_k(k, R);
    }
    PairEnergyCOM(){}

    double all2all(Particles& particles){
        double e = 0.0;

        //printf("all2all geo: %lf %lf %lf\n", this->geo->dh[0], this->geo->dh[1], this->geo->dh[2]);
        #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = i + 1; j < particles.tot; j++){
                //printf("distance %lf\n", this->geo->distance(particles[i]->pos, particles[j]->pos));
                //printf("Indices: %u, %u\n", i, j);
                e += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->com, particles[j]->com));
            }  
        }
        //printf("Real energy: %.15lf\n", e);
        return e;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = 0.0;
        
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            e += i2i(p->q, particles[i]->q, this->geo->distance(p->com, particles[i]->com));
        }

        return e;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        //printf("i2all geo: %lf %lf %lf\n", this->geo->dh[0], this->geo->dh[1], this->geo->dh[2]);
        double e = 0.0;
        // Need to fix this, not a nice solution (when volume move).........
        if(p.size() == particles.tot){
            e = all2all(particles);
        }
        else{
            #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
            for(int i = 0; i < p.size(); i++){
                //do instead i2all(s, particles);
                e += i2all(particles.particles[p[i]], particles);
            }

            for(int i = 0; i < p.size(); i++){
                for(int j = i + 1; j < p.size(); j++){
                    e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->com, particles[p[j]]->com));
                }
            }
        }

        return e;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){
        //printf("i2all geo: %lf %lf %lf\n", this->geo->dh[0], this->geo->dh[1], this->geo->dh[2]);
        double e = 0.0;
        // Need to fix this, not a nice solution.........
        if(p.size() == particles.tot){
            e = all2all(particles);
        }
        else{
            #pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
            for(int i = 0; i < p.size(); i++){
                //do instead i2all(s, particles);
                e += i2all(particles.particles[p[i]], particles);
            }

            for(int i = 0; i < p.size(); i++){
                for(int j = i + 1; j < p.size(); j++){
                    e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->com, particles[p[j]]->com));
                }
            }
        }

        return e;
    }

    inline double i2i(double& q1, double& q2, double&& dist){
        if(dist <= this->cutoff){
            return energy_func(q1, q2, dist);
        }
        else{
            return 0.0;
        }
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void initialize(Particles& particles){}
    void update(double x, double y, double z){}

    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for(unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            force += energy_func.force(p->q, particles[i]->q, this->geo->displacement(p->com, particles[i]->com));
        }
        return force;
    }

    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force;
        return force(p1->com, p2->com, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d disp = this->geo->displacement(pos1, pos2);
        if(disp.norm() <= this->cutoff){
            return energy_func.force(q1, q2, disp);
        }
        else{
            return Eigen::Vector3d::Zero();
        }

    }
};











template <typename E>
class ChargeWell : public EnergyBase{
    private:

    E energy_func;  //energy functor

    public:
    ChargeWell(double k, double R){
        energy_func.set_k(k, R);
    }

    inline double i2i(double& r, double&& dist){
        if(dist <= this->cutoff){
            return energy_func(r, dist);
        }
        else{
            return 0.0;
        }
    }

    double all2all(Particles& particles){
        double e = 0.0;
        for(unsigned int i = 0; i < particles.tot; i++){
            e += this->i2i(particles[i]->r, this->geo->distance(particles[i]->pos, particles[i]->com));
        }
        return e;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = this->i2i(p->r, this->geo->distance(p->pos, p->com));
        return e;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e;
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void initialize(Particles& particles){}
    void update(double x, double y, double z){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force;
        return force(p1->pos, p2->pos, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d force;
        force = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));
        return force;
    }
};





template <typename E>
class ExternalEnergy : public EnergyBase{
    private:

    E energy_func;  //energy functor

    public:
    ExternalEnergy(double x, double y, double z, double k){
        energy_func.set_bounds(x, y, z, k);
    }

    inline double i2i(double& q, Eigen::Vector3d& disp){
        return energy_func(q, disp);
    }

    double all2all(Particles& particles){
        double e = 0.0;

        for(unsigned int i = 0; i < particles.tot; i++){
            e += this->i2i(particles[i]->q, particles[i]->pos);
        }

        return e * constants::lB;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){

        double e = this->i2i(p->q, p->pos);;

        return e;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }

        return e * constants::lB;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e * constants::lB;
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void initialize(Particles& particles){}
    void update(double x, double y, double z){}

    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){
        return energy_func.force(p->q, p->pos) * constants::lB;
    }

    Eigen::Vector3d force(std::shared_ptr<Particle> p, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force * constants::lB;
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        return Eigen::Vector3d::Zero();
    }
};








template <typename E>
class PairEnergyWithRep : public EnergyBase{

    private:

    E energy_func;  //energy functor
    int rep;
    public:

    PairEnergyWithRep(int rep){
        this->rep = rep;
    }


    double all2all(Particles& particles){
        double e = 0.0;
        Eigen::Vector3d disp;

        //#pragma omp parallel for reduction(+:e) schedule(guided, 100) if(particles.tot >= 500)
        for(int l = -this->rep; l <= this->rep; l++){
            for(int m = -this->rep; m <= this->rep; m++){
                for(int n = -this->rep; n <= this->rep; n++){
                    for(unsigned int i = 0; i < particles.tot; i++){
                        for(unsigned int j = 0; j < particles.tot; j++){
                            if(l == 0 && m == 0 && n == 0 && i == j) continue;
                            disp << particles[j]->pos[0] + l * geo->d[0], particles[j]->pos[1] +  m * geo->d[1], particles[j]->pos[2] +  n * geo->d[2];
                            e += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, disp));
                        }  
                    }
                }
            }
        }

        e *= 0.5;
        printf("Real energy: %.15lf\n", e);
        return e * constants::lB;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = 0.0;

        Eigen::Vector3d disp;

        //#pragma omp parallel for reduction(+:e) schedule(guided, 100) if(particles.tot >= 500)
        for(int l = -this->rep; l <= this->rep; l++){
            for(int m = -this->rep; m <= this->rep; m++){
                for(int n = -this->rep; n <= this->rep; n++){
                    for (unsigned int i = 0; i < particles.tot; i++){
                        if (l == 0 && m == 0 && n == 0 && p->index == particles[i]->index) continue;
                        disp << particles[i]->pos[0] + (double)l * geo->d[0], particles[i]->pos[1] +  (double)m * geo->d[1], particles[i]->pos[2] +  (double)n * geo->d[2];
                        if(p->index == particles[i]->index){
                            e += 0.5 * i2i(p->q, particles[i]->q, this->geo->distance(p->pos, disp));
                        }
                        else{
                            e += i2i(p->q, particles[i]->q, this->geo->distance(p->pos, disp));
                        }    
                    }
                }
            }
        }

        return e;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            //do instead i2all(s, particles);
            e += i2all(particles.particles[s], particles);
        }

        for(int i = 0; i < p.size(); i++){
           for(int j = i + 1; j < p.size(); j++){
               e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
           }
        }

        return e * constants::lB;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            //do instead i2all(s, particles);
            e += i2all(particles.particles[s], particles);
        }

        for(int i = 0; i < p.size(); i++){
           for(int j = i + 1; j < p.size(); j++){
               e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
           }
        }

        return e * constants::lB;
    }

    inline double i2i(double& q1, double& q2, double&& dist){
        if(dist <= this->cutoff){
            return energy_func(q1, q2, dist);
        }
        else{
            return 0.0;
        }
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void update(double x, double y, double z){}
    void initialize(Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force * constants::lB;
        return force(p1->pos, p2->pos, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d force;
        force = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));
        return force * constants::lB;
    }
};











template <typename E>
class ExtEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:
    ExtEnergy(double x, double y, double z){
        energy_func.set_box(x, y, z);
    }

    double i2all(std::shared_ptr<Particle> p, Particles& particles){ 
        return energy_func();
    }

    double i2i(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){ return 0.0; }

    double all2all(Particles& particles){
        //printf("Rec energy: %.15lf\n", energy_func());
        return energy_func() * constants::lB;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        return energy_func() * constants::lB;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){
        return energy_func() * constants::lB;
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){
        energy_func.update(_old, _new);
    }

    void update(double x, double y, double z){
        energy_func.set_box(x, y, z);
    }

    void initialize(Particles& particles){
        energy_func.initialize(particles);
    }

    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for(unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            force += energy_func.force(p->q, particles[i]->q, this->geo->displacement(p->pos, particles[i]->pos));
        }
        return force * constants::lB;
    }

    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        Eigen::Vector3d force;
        force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        return force * constants::lB;
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d force;
        force = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));
        return force * constants::lB;
    }
};




template <typename E>
class ImgEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double CC = 0.0, CpC = 0.0, self = 0.0;
        Eigen::Vector3d temp;

        // CC
        #pragma omp parallel for reduction(+:CC) schedule(guided, 500) if(particles.tot >= 3000) 
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;

            CC += i2i(p->q, particles[i]->q, this->geo->distance(p->pos, particles[i]->pos));
        }


        //  C'C
        temp = p->pos;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        #pragma omp parallel for reduction(+:CpC) schedule(guided, 500) if(particles.tot >= 3000) 
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;

            CpC += i2i(-p->q, particles[i]->q, this->geo->distance(temp, particles[i]->pos));
        }
        // => CC == C'C' and C'C == CC'

        // Self term  

        self = i2i(p->q, -p->q, this->geo->distance(p->pos, temp));
        return CC + CpC + 0.5 * self;
    }


    inline double i2i(double q1, double q2, double&& dist){
        if(dist <= this->cutoff){
            return energy_func(q1, q2, dist);
        }
        else{
            return 0.0;
        }
    }


    double all2all(Particles& particles){
        double CC = 0.0, CpC = 0.0;
        Eigen::Vector3d temp;
        Eigen::Vector3d temp2;

        // CC
        #pragma omp parallel for schedule(guided, 200) reduction(+:CC) if(particles.tot >= 1000)
        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = i + 1; j < particles.tot; j++){
                //CC += energy_func(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
                CC += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
            } 
        }

        //C'C
        #pragma omp parallel for schedule(dynamic, 200) reduction(+:CpC) private(temp) if(particles.tot >= 1000)
        for(unsigned int i = 0; i < particles.tot; i++){
            temp = particles[i]->pos;
            temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 

            for(unsigned int j = 0; j < particles.tot; j++){
                //CpC += energy_func(-particles[i]->q, particles[j]->q, this->geo->distance(temp, particles[j]->pos));
                CpC += i2i(-particles[i]->q, particles[j]->q, this->geo->distance(temp, particles[j]->pos));
            } 
        }

        return (CC + 0.5 * CpC) * constants::lB;
    }


    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles[s], particles);
        }

        for(int i = 0; i < p.size(); i++){
           for(int j = i + 1; j < p.size(); j++){
               e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
           }
        }

        return e * constants::lB;
    }


    double operator()(std::vector< unsigned int >& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles[s], particles);
        }

        for(int i = 0; i < p.size(); i++){
           for(int j = i + 1; j < p.size(); j++){
               e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
           }
        }

        return e * constants::lB;
    }


    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){
        UNUSED(_old);
        UNUSED(_new);
    }

    void update(double x, double y, double z){}

    void initialize(Particles& particles){
        UNUSED(particles);
    }


    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for(unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            force += energy_func.force(p->q, particles[i]->q, this->geo->displacement(p->pos, particles[i]->pos));
        }
        return force * constants::lB;

        Eigen::Vector3d CC = Eigen::Vector3d::Zero();
        Eigen::Vector3d CpC = Eigen::Vector3d::Zero();;
        Eigen::Vector3d self = Eigen::Vector3d::Zero();;
        Eigen::Vector3d temp;

        // CC
        #pragma omp parallel for reduction(+:CC) schedule(guided, 500) if(particles.tot >= 3000) 
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;

            CC += energy_func.force(p->q, particles[i]->q, this->geo->displacement(p->pos, particles[i]->pos));
        }


        //  C'C
        temp = p->pos;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        #pragma omp parallel for reduction(+:CpC) schedule(guided, 500) if(particles.tot >= 3000) 
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;

            CpC += energy_func.force(-p->q, particles[i]->q, this->geo->displacement(temp, particles[i]->pos));
        }
        // => CC == C'C' and C'C == CC'

        // Self term  

        self = energy_func.force(p->q, -p->q, this->geo->displacement(p->pos, temp));
        return (CC - CpC + 0.5 * self) * constants::lB;
    }


    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        Eigen::Vector3d CC;
        Eigen::Vector3d CpC;
        Eigen::Vector3d self;
        Eigen::Vector3d temp;

        // CC
        CC = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));

        //  C'C
        temp = p2->pos;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        CpC = energy_func.force(p1->q, -p2->q, this->geo->displacement(p1->pos, temp));
        
        // Self term  
        temp = p1->pos;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        self = energy_func.force(p1->q, -p1->q, this->geo->displacement(p1->pos, temp));

        //return 0.5 * (CC + CpC + 0.5 * self);
        return CC * constants::lB;
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d CC;
        Eigen::Vector3d CpC;
        Eigen::Vector3d self;
        Eigen::Vector3d temp;

        // CC
        CC = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));

        //  C'C
        temp = pos2;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        CpC = energy_func.force(q1, -q2, this->geo->displacement(pos1, temp));
        
        // Self term  
        temp = pos1;
        temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
        self = energy_func.force(q1, -q2, this->geo->displacement(pos1, temp));

        //return 0.5 * (CC + CpC + 0.5 * self);
        return CC * constants::lB;
    }
};






template <typename E>
class MIHalfwald : public EnergyBase{

    private:

    E energy_func;  //energy functor
    int kMax;
    double eps;

    public:

    MIHalfwald(int kMax, double eps) : kMax(kMax), eps(eps){
        printf("\tNumber of replicas (on each side of original cell): %i\n", this->kMax);
        printf("\teps factor: %lf\n", this->eps);
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double CC = 0.0, CpC = 0.0;
        //double self = 0.0;
        Eigen::Vector3d temp;

        // CC
        for(int k = -this->kMax; k <= this->kMax; k++){
            #pragma omp parallel for schedule(dynamic, 250) reduction(+:CC) private(temp) if(particles.tot >= 1000)
            for (unsigned int i = 0; i < particles.tot; i++){
                if (p->index == particles[i]->index && k == 0) continue;
                double tmpE = 0.0;

                temp = particles[i]->pos;
                temp[2] += k * 2.0 * this->geo->_d[2]; 
                tmpE = i2i(p->q, particles[i]->q, this->geo->distance(p->pos, temp)) * std::pow(this->eps, 2.0 * std::fabs(k)); 

                if(p->index == i){
                    tmpE *= 0.5; 
                }
                
                CC += tmpE;
            }
        }

        //  CC'
        for(int k = -this->kMax; k <= this->kMax; k++){
            #pragma omp parallel for schedule(dynamic, 250) reduction(+:CpC) private(temp) if(particles.tot >= 1000)
            for (unsigned int i = 0; i < particles.tot; i++){
                double tmpE = 0.0;
                temp = particles[i]->pos;
                //temp[2] = math::sgn(temp[2]) * this->geo->dh[2] - temp[2]; 
                temp[2] = math::sgn(temp[2]) * this->geo->_d[2] - temp[2] + k * 2.0 * this->geo->_d[2]; 

                tmpE = i2i(p->q, -particles[i]->q, this->geo->distance(p->pos, temp)) * std::pow(this->eps, 2.0 * std::fabs(k) + 1.0);
                if(p->index == i){
                    tmpE *= 0.5; 
                }

                CpC += tmpE;
            }
        }
        // => CC == C'C' and C'C == CC'

        return CC + CpC;
    }


    inline double i2i(double q1, double q2, double&& dist){
        if(dist <= this->cutoff){
            return energy_func(q1, q2, dist);
        }
        else{
            return 0.0;
        }
    }


    double all2all(Particles& particles){
        double CC = 0.0, CpC = 0.0;
        Eigen::Vector3d temp;

        // CC box-box
        for(int k = -this->kMax; k <= this->kMax; k++){
            #pragma omp parallel for schedule(guided, 200) reduction(+:CC) private(temp) if(particles.tot >= 1000)
            for(unsigned int i = 0; i < particles.tot; i++){
                for(unsigned int j = 0; j < particles.tot; j++){
                    if(k == 0 && i == j) continue;

                    double tmpE = 0.0;
                    temp = particles[j]->pos;
                    temp[2] += k * 2.0 * this->geo->_d[2]; 

                    tmpE = i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, temp)) * std::pow(this->eps, 2.0 * std::fabs(k));

                    CC += tmpE;
                } 
            }
        }


        //CC'
        for(int k = -this->kMax; k <= this->kMax; k++){
            #pragma omp parallel for schedule(dynamic, 200) reduction(+:CpC) private(temp) if(particles.tot >= 1000)
            for(unsigned int i = 0; i < particles.tot; i++){
                for(unsigned int j = 0; j < particles.tot; j++){
                    double tmpE = 0.0;
                    temp = particles[j]->pos;
                    temp[2] = math::sgn(temp[2]) * this->geo->_d[2] - temp[2] + k * 2.0 * this->geo->_d[2];
                    tmpE = i2i(particles[i]->q, -particles[j]->q, this->geo->distance(particles[i]->pos, temp)) * std::pow(this->eps, 2.0 * std::fabs(k) + 1.0);

                    CpC += tmpE;
                } 
            }
        }

        return (0.5 * CC + 0.5 * CpC) * constants::lB;
    }


    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles[s], particles);
        }

        for(int i = 0; i < p.size(); i++){
           for(int j = i + 1; j < p.size(); j++){
               e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
           }
        }

        return e * constants::lB;
    }


    double operator()(std::vector< unsigned int >& p, Particles& particles){
        double e = 0.0;
        for(auto s : p){
            e += i2all(particles[s], particles);
        }

        for(int i = 0; i < p.size(); i++){
           for(int j = i + 1; j < p.size(); j++){
               e -= i2i(particles[p[i]]->q, particles[p[j]]->q, this->geo->distance(particles[p[i]]->pos, particles[p[j]]->pos));
           }
        }

        return e * constants::lB;
    }


    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){
        UNUSED(_old);
        UNUSED(_new);
    }

    void update(double x, double y, double z){}

    void initialize(Particles& particles){
        UNUSED(particles);
    }
    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force * constants::lB;
        return force(p1->pos, p2->pos, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d force;
        force = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));
        return force * constants::lB;
    }
};




template <typename E>
class Ellipsoid : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    Ellipsoid(std::vector<double> a, std::vector<double> b, std::vector<double >c){
        energy_func.load(a, b, c);
    }

    double all2all(Particles& particles){
        double e = 0.0;


        //#pragma omp parallel for reduction(+:e) schedule(guided, 100) if(particles.tot >= 500)
        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = i + 1; j < particles.tot; j++){
                Eigen::Vector3d disp = particles[i]->pos - particles[j]->pos;
                e += i2i(particles[i]->q, particles[j]->q, disp);
            } 
        }
        printf("Real energy: %.15lf\n", e);
        return e * constants::lB;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = 0.0;

        //#pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            Eigen::Vector3d disp = p->pos - particles[i]->pos;
            e += i2i(p->q, particles[i]->q, disp);
        }

        return e;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }

        return e * constants::lB;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){

        double e = 0.0;
        for(auto s : p){
            e += i2all(particles.particles[s], particles);
        }
        return e * constants::lB;
    }

    inline double i2i(double q1, double q2, Eigen::Vector3d disp){
        double a = 50.0, b = 50.0, c = 25.0;
        if(disp[0] * disp[0] / (a * a) + disp[1] * disp[1] / (b * b) + disp[2] * disp[2] / (c * c) <= 1.0){
            return energy_func(q1, q2, disp);
        }
        else{
            return 0.0;
        }
    }

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void update(double x, double y, double z){}
    void initialize(Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force * constants::lB;
        return force(p1->pos, p2->pos, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d force;
        force = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));
        return force * constants::lB;
    }
};













template <typename E>
class Energy2D: public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:
    Energy2D(double x, double y, double z){
        energy_func.set_box(x, y, z);
    }

    double all2all(Particles& particles){
        double reciprocal = 0.0, real = 0.0, self = 0.0, charge = 0.0, gE = 0.0, gEt = 0.0;
        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = 0; j < particles.tot; j++){
                if(i == j) continue;
                real += energy_func(particles.particles[i], particles.particles[j], this->geo->distance(particles[i]->pos, particles[j]->pos));
            }  
            charge += particles[i]->q;
        }

        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = 0; j < particles.tot; j++){
                //Include i=j term
                reciprocal += energy_func.rec(particles[i], particles[j], this->geo->displacement(particles[i]->pos, particles[j]->pos));
                gEt = energy_func.gE(particles[i]->q, particles[j]->q, this->geo->displacement(particles[i]->pos, particles[j]->pos));
                gE += gEt;
            }  
        }

        self = energy_func.get_self();
        real *= 0.5;
        //reciprocal *= 0.5;
        //gE *= 0.5;
        printf("real %lf\n", real * constants::lB);
        printf("self %lf\n", self * constants::lB);
        printf("reciprocal %lf\n", reciprocal * constants::lB);
        printf("gE %lf\n", gE * constants::lB);
        double energy = (real + reciprocal - gE - self);
        //printf("Real energy: %.15lf\n", e);
        return energy * constants::lB;
    }

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double e = all2all(particles) / constants::lB;
        return e * constants::lB;
    }

    double operator()(std::vector< unsigned int >&& p, Particles& particles){
        double e = all2all(particles) / constants::lB;
        return e * constants::lB;
    }

    double operator()(std::vector< unsigned int >& p, Particles& particles){
        double e = all2all(particles) / constants::lB;
        return e * constants::lB;
    }

    inline double i2i(double& q1, double& q2, double&& dist){}

    void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new){}
    void initialize(Particles& particles){
        energy_func.initialize(particles);
    }
    void update(double x, double y, double z){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p, Particles& particles){}
    Eigen::Vector3d force(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
        //Eigen::Vector3d force;
        //force = energy_func.force(p1->q, p2->q, this->geo->displacement(p1->pos, p2->pos));
        //return force * constants::lB;
        return force(p1->pos, p2->pos, p1->q, p2->q);
    }

    Eigen::Vector3d force(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double q1, double q2){
        Eigen::Vector3d force;
        force = energy_func.force(q1, q2, this->geo->displacement(pos1, pos2));
        return force * constants::lB;
    }
};