#pragma once

#include "particle.h"
#include "particles.h"
#include "geometry.h"


class EnergyBase{

    protected:

    Geometry *geo;
    double cutoff;

    public:

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
    virtual double operator()(std::vector< unsigned int >& p, Particles& particles) = 0;
    virtual void update(std::vector< std::shared_ptr<Particle> >&& _old, std::vector< std::shared_ptr<Particle> >&& _new) = 0;
    virtual void initialize(Particles& particles) = 0;
};


template <typename E>
class PairEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:

    double all2all(Particles& particles){
        double e = 0.0;


        //#pragma omp parallel for reduction(+:e) schedule(guided, 100) if(particles.tot >= 500)
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

        //#pragma omp parallel for reduction(+:e) schedule(dynamic, 100) if(particles.tot >= 500)
        for (unsigned int i = 0; i < particles.tot; i++){
            if (p->index == particles[i]->index) continue;
            e += i2i(p->q, particles[i]->q, this->geo->distance(p->pos, particles[i]->pos));
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
    void initialize(Particles& particles){}
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
    void initialize(Particles& particles){}
};











template <typename E>
class ExtEnergy : public EnergyBase{

    private:

    E energy_func;  //energy functor

    public:
    ExtEnergy(double x, double y, double z){
        energy_func.set_box(x, y, z);
    }

    double i2all(std::shared_ptr<Particle> p, Particles& particles){ return 0.0; }
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
    void initialize(Particles& particles){
        energy_func.initialize(particles);
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


    inline double i2i(const double q1, const double q2, const double&& dist){
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


    void initialize(Particles& particles){
        UNUSED(particles);
    }


};






template <typename E>
class MIHalfwald : public EnergyBase{

    private:

    E energy_func;  //energy functor
    int kMax;

    public:

    MIHalfwald(int kMax) : kMax(kMax){}

    inline double i2all(std::shared_ptr<Particle> p, Particles& particles){
        double CC = 0.0, CpC = 0.0, self = 0.0;
        Eigen::Vector3d temp;

        // CC
        for(int k = -this->kMax; k <= this->kMax; k++){
            //#pragma omp parallel for reduction(+:CC) schedule(guided, 500) if(particles.tot >= 3000) 
            for (unsigned int i = 0; i < particles.tot; i++){
                if (p->index == particles[i]->index && k == 0) continue;

                temp = particles[i]->pos;
                temp[2] += k * this->geo->d[2]; 

                if (p->index == particles[i]->index){
                    CC += 0.5 * i2i(p->q, particles[i]->q, this->geo->distance(p->pos, temp));
                }
                else{
                    CC += i2i(p->q, particles[i]->q, this->geo->distance(p->pos, temp));    
                }
            }
        }

        //  CC'
        for(int k = -this->kMax; k <= this->kMax; k++){
            //#pragma omp parallel for reduction(+:CpC) schedule(guided, 500) if(particles.tot >= 3000) 
            for (unsigned int i = 0; i < particles.tot; i++){
                temp = particles[i]->pos;
                temp[2] = math::sgn(temp[2]) * this->geo->_d[2] - temp[2] + k * this->geo->d[2]; 
                if (p->index == particles[i]->index){
                    CpC += 0.5 * i2i(p->q, -particles[i]->q, this->geo->distance(p->pos, temp));
                }
                else{
                    CpC += i2i(p->q, -particles[i]->q, this->geo->distance(p->pos, temp));
                }
            }
        }
        // => CC == C'C' and C'C == CC'

        return CC + CpC;
    }


    inline double i2i(const double q1, const double q2, const double&& dist){
        if(dist <= this->cutoff){
            return energy_func(q1, q2, dist);
        }
        else{
            return 0.0;
        }
    }


    double all2all(Particles& particles){
        double CC = 0.0, CpC = 0.0, CCr = 0.0;
        Eigen::Vector3d temp;

        // CC box-box
        //#pragma omp parallel for schedule(guided, 200) reduction(+:CC) if(particles.tot >= 1000)
        for(unsigned int i = 0; i < particles.tot; i++){
            for(unsigned int j = i + 1; j < particles.tot; j++){
                //CC += energy_func(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
                CC += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
            } 
        }

        // CC box-replicates
        for(int k = -this->kMax; k <= this->kMax; k++){
            if(k == 0) continue;
            //#pragma omp parallel for schedule(guided, 200) reduction(+:CC) if(particles.tot >= 1000)
            for(unsigned int i = 0; i < particles.tot; i++){
                for(unsigned int j = 0; j < particles.tot; j++){
                    temp = particles[j]->pos;
                    temp[2] += k * this->geo->d[2]; 
                    //CC += energy_func(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, particles[j]->pos));
                    if(i == j){
                        CCr += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, temp));
                    }
                    else{
                        CCr += i2i(particles[i]->q, particles[j]->q, this->geo->distance(particles[i]->pos, temp));
                    }
                } 
            }
        }


        //CC'
        for(int k = -this->kMax; k <= this->kMax; k++){
            //#pragma omp parallel for schedule(dynamic, 200) reduction(+:CpC) private(temp) if(particles.tot >= 1000)
            for(unsigned int i = 0; i < particles.tot; i++){
                for(unsigned int j = 0; j < particles.tot; j++){
                    temp = particles[j]->pos;
                    temp[2] = math::sgn(temp[2]) * this->geo->_d[2] - temp[2] + k * this->geo->d[2]; 
                    //CpC += energy_func(-particles[i]->q, particles[j]->q, this->geo->distance(temp, particles[j]->pos));
                    if(i == j){
                        CpC += i2i(particles[i]->q, -particles[j]->q, this->geo->distance(particles[i]->pos, temp));
                    }

                    else{
                        CpC += i2i(particles[i]->q, -particles[j]->q, this->geo->distance(particles[i]->pos, temp));
                    }
                } 
            }
        }

        return (CC + 0.5 * CpC + 0.5 * CCr) * constants::lB;
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


    void initialize(Particles& particles){
        UNUSED(particles);
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
    void initialize(Particles& particles){}
};
