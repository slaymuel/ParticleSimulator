#include "particle.h"
#include <vector>
//#include "geometry.h"
#include "Faddeeva.h"

/*
#pragma omp declare reduction(vec_double_plus : std::vector<std::complex<double>> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<std::complex<double>>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
*/


class Coulomb{
    public:

    inline double operator()(const double& q1, const double& q2, const double& dist){
        return q1 * q2 / dist;
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }
};


class LJ{
    private:
    const double s = 5.0;
    double s6 = s*s*s*s*s*s;
    double s12 = s6*s6;
    double k = 1.0;

    public:
    void set_k(double k, double R){
        this->k = k;
    }

    inline double operator()(const double& q1, const double& q2, double& dist){
        double d6 = dist*dist*dist*dist*dist*dist;

        return 4.0 * k * (s12 / (d6 * d6) - s6 / d6);
    }

    inline Eigen::Vector3d force(double& q1, double& q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        double n = disp.norm();
        double d6 = n*n*n*n*n*n;
        force = 24.0 * k * (2.0 * s12 / (d6 * d6 * n) - s6 / (d6 * n)) * disp.normalized();
        
        return force;
    }
};

class LJST{
    private:
    const double s = 5.0;

    double k = 1.0;
    double shift;
    double s6 = s*s*s*s*s*s;
    double s12 = s6*s6;

    public:
    void set_k(double k, double R){
        this->shift = 0.0;
        this->k = k;
        this->shift = (*this)(1.0, 1.0, R);
        printf("\tShift: %lf", this->shift);
    }

    inline double operator()(const double& q1, const double& q2, double& dist){
        double d6 = dist*dist*dist*dist*dist*dist;

        return 4.0 * k * (s12 / (d6 * d6) - s6 / d6) - this->shift;
    }

    inline Eigen::Vector3d force(double& q1, double& q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        double n = disp.norm();
        double d6 = n*n*n*n*n*n;
        force = 24.0 * k * (2.0 * s12 / (d6 * d6 * n) - s6 / (d6 * n)) * disp.normalized();
        
        return force;
    }
};


class LJRep{
    private:
    const double s = 4.0;
    double s6 = s*s*s*s*s*s;
    double s12 = s6*s6;
    double k;

    public:

    LJRep(){
        this->k = 1.0; //constants::KB * constants::T;
    }
    inline double operator()(const double& q1, const double& q2, double& dist){
        double d6 = dist*dist*dist*dist*dist*dist;
        
        return this->k * s12 / (d6 * d6);
    }

    inline Eigen::Vector3d force(double& q1, double& q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        double n = disp.norm();
        double d6 = n*n*n*n*n*n;
        force = -12.0 * this->k * s12 * s / (d6 * d6 * n) * disp.normalized();
        
        return force;
    }
};



class LJWallRep{
    private:
    double k;
    double zWall;

    public:
    void set_bounds(double x, double y, double z, double k){
        this->zWall = 0.5 * z;
        this->k = k;
    }

    inline double operator()(const double& q, const Eigen::Vector3d& p){
        double energy;

        //Left wall
        double dist = p[2] + zWall;
        double d6 = dist*dist*dist*dist*dist*dist;
        energy =  1.0 / (d6 * d6);
        //Right wall
        dist = 2.0 * zWall - dist;
        d6 = dist*dist*dist*dist*dist*dist;
        energy += 1.0 / (d6 * d6);

        return this->k * energy;
    }

    inline Eigen::Vector3d force(const double& q, const Eigen::Vector3d& p){
        Eigen::Vector3d force;

        //Left wall
        double dist = p[2] + zWall;
        double d6 = dist*dist*dist*dist*dist*dist;
        force =  1.0 / (d6 * d6 * dist) * Eigen::Vector3d::UnitZ();
        //Right wall
        dist = 2.0 * zWall - dist;
        d6 = dist*dist*dist*dist*dist*dist;
        force -= 1.0 / (d6 * d6 * dist) * Eigen::Vector3d::UnitZ();

        return -12.0 * this->k * force;
    }
};


class ExpWallRep{
    private:
    float k;
    float zWall;

    public:
    void set_bounds(double x, double y, double z, double k){
        this->zWall = 0.5 * z;
        this->k = k;
    }

    inline double operator()(const double& q, const Eigen::Vector3d& p){
        double energy;

        //Left wall
        double dist = p[2] + zWall;
        energy = std::exp(-dist);
        //Right wall
        dist = 2.0 * zWall - dist;
        energy += std::exp(-dist);

        return this->k * energy;
    }

    inline Eigen::Vector3d force(const double& q, const Eigen::Vector3d& p){
        Eigen::Vector3d force;

        //Left wall
        double dist = p[2] + zWall;
        double d6 = dist*dist*dist*dist*dist*dist;
        force =  1.0 / (d6 * d6 * dist) * Eigen::Vector3d::UnitZ();
        //Right wall
        dist = 2.0 * zWall - dist;
        d6 = dist*dist*dist*dist*dist*dist;
        force -= 1.0 / (d6 * d6 * dist) * Eigen::Vector3d::UnitZ();

        return -12.0 * this->k * force;
    }
};


class Harmonic{
    private:
    double k;

    public:
    void set_k(double k, double R){
        this->k = k;
        printf("\tForce constant is (k): %lf\n", this->k);
    }

    inline double operator()(const double& R, const double& dist){
        //printf("k: %lf dist: %lf\n", this->k, dist);
        return this->k * dist * dist;
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }
};



class FENE{
    private:
    double k, Rsq;

    public:
    void set_k(double k, double R){
        this->k = k;
        this->Rsq = R * R;
        printf("\tForce constant is (k): %lf, R is: %lf\n", this->k, R);
    }

    inline double operator()(const double& R, const double& dist){
        if(dist * dist > this->Rsq) return 1e30;
        return -0.5 * this->k * this->Rsq * std::log(1.0 - dist * dist / this->Rsq);
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }
};




class Jan{
    private:
    double k, R;

    public:
    void set_k(double k, double R){
        this->k = k;
        this->R = R;
        printf("\tForce constant is (k): %lf, R is: %lf\n", this->k, R);
    }

    inline double operator()(const double& R, const double& dist){
        double sq = (this->R - dist) * (this->R - dist);
        return this->k * dist * dist * (this->R * this->R / sq);
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }
};




class Sture{
    private:
    double k, Rsq;

    public:
    void set_k(double k, double R){
        this->k = k;
        this->Rsq = R * R;
        printf("\tForce constant is (k): %lf, R is: %lf\n", this->k, R);
    }

    inline double operator()(const double& R, const double& dist){
        if(dist * dist >= this->Rsq) return 10e250;//std::numeric_limits<double>::infinity();

        return this->k * this->Rsq * (this->Rsq / (this->Rsq - dist * dist) - 1.0);
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }
};




namespace Fanourgakis{
    double R;

    class SP2{
        private:

        public:

        inline double operator()(const double& q1, const double& q2, const double& dist){
            return q1 * q2 * (1.0 - 2.0 * dist / R + 2.0 * (dist / R) * (dist / R) * (dist / R) - (dist / R) * (dist / R) * 
                                                                                                  (dist / R) * (dist / R)) / dist;
        }

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };



    class SP2Self{
        private:

        double selfTerm = 0.0;

        public:

        void set_box(double x, double y, double z){}



        void initialize(Particles &particles){
            this->selfTerm = 0.0;

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
            printf("Self term: %lf\n", this->selfTerm);
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q;
                }
            }

            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q;
                }
            }
        }


        inline double operator()(){
            return -1.0 / R * this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };


    class SP3{
        private:

        public:

        inline double operator()(const double& q1, const double& q2, const double& dist){
            return q1 * q2 * (1.0 - 7.0 / 4.0 * dist / R + 21.0 / 4.0 * std::pow((dist / R), 5.0) - 
                   7.0 * std::pow((dist / R), 6.0) + 5.0 / 2.0 * std::pow((dist / R), 7.0)) / dist;
        }
        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };


    class SP3Self{
        private:

        double selfTerm = 0.0;

        public:

        void set_box(double x, double y, double z){}



        void initialize(Particles &particles){
            this->selfTerm = 0.0;
            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q;
                }
            }

            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q;
                }
            }
        }


        inline double operator()(){
            return -7.0 / (8.0 * R) * this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };
}





class Spline{
    std::vector<double> knots;

    public:
    Spline(std::vector<double> k) : knots(k) {}

    inline double spline(double t, int i, int k){
        double value = 0.0;

        if(k == 0){ //End condition
            if(knots[i] <= t && t < knots[i + 1]){
                return 1.0;
            }
            else{
                return 0.0;
            }
        }
        else{
            if(knots[i + k] - knots[i] != 0.0 && knots[i + k + 1] -  knots[i + 1] != 0.0){
                value += (t - knots[i]) / (knots[i + k] - knots[i]) * this->spline(t, i, k - 1) + (knots[i + k + 1] - t) / (knots[i + k + 1] -  knots[i + 1]) * this->spline(t, i + 1, k - 1);
            }
        }

        return value;
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }
};




class BSpline2D{


    std::vector<double> aKnots;
    std::vector<double> bKnots;
    std::vector<double> controlPoints;
    Spline *aSpline;
    Spline *bSpline;
    int degree;

    public:
    void load(std::vector<double> aKnots, std::vector<double> bKnots, std::vector<double >controlPoints){
        aSpline = new Spline(aKnots);
        bSpline = new Spline(bKnots);
        this->controlPoints = controlPoints;
        degree = 3;
        printf("\tSpline loaded into potential, elements: %lu, %lu, %lu\n", aKnots.size(), bKnots.size(), controlPoints.size());
    }


    inline double operator()(const double& q1, const double& q2, Eigen::Vector3d& disp){
        return q1 * q2 * get_energy(disp);
    }


    inline double get_energy(Eigen::Vector3d& disp){
        double energy = 0.0;
        int max = std::sqrt(controlPoints.size());
        double p = std::sqrt(disp[0] * disp[0] + disp[1] * disp[1]);

        for(int i = 0; i < max; i++){
            for(int j = 0; j < max; j++){
                energy += controlPoints[i * max + j] * aSpline->spline(p, i, degree) * bSpline->spline(disp[2], j, degree);
            }
        }

        return energy;
    }

    inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
        Eigen::Vector3d force;
        return force;
    }

};



















namespace EwaldLike{
    double alpha = 0.0;
    double kMax = 0.0;
    std::vector<int> kM;
    double R;
    double eta;
    bool spherical;

    void set_km(std::vector<int> v){
        kM = v;
    }




    class Short{

        public:

        inline double operator()(const double& q1, const double& q2, const double& dist){

            //math::sgn(p2->pos[2]) * d[2] - p2->pos[2];   //Mirror of p2
            double energy = q1 * q2 / dist;
            double real = math::erfc_x(dist * alpha) * energy;

            //printf("Real %.15lf\n", real);
            return real;    //tinfoil
        }

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d real;
            real << 0.0, 0.0, 0.0;
            double dist = disp.norm();

            //Real force
            real = disp / (dist * dist);
            real *= 2.0 * std::sqrt(alpha / constants::PI) * std::exp(-alpha * dist*dist) + 1.0 / dist * math::erfc_x(std::sqrt(alpha) * dist);
            /*printf("\nReal force: \n");
            std::cout << q1 * q2 * real << std::endl;
            printf("\n");*/
            return q1 * q2 * real;
        }
    };





    class ShortTruncated{
        private:

        public:
        inline double operator()(const double& q1, const double& q2, const double& dist){
            double energy = 0.0;
            double q = dist / R;

            if(dist >= R){
                return 0.0;
            }
            else if(dist < 1e-6){
                printf("Distance is 0\n");
                exit(0);
            }
            else{
                //energy = math::erfc_x(R * std::sqrt(2.0) * q / (2.0 * alpha)) - math::erfc_x(R * std::sqrt(2) / (2.0 * alpha)) - (1.0 - q) * R * std::sqrt(2.0) * std::exp(-R*R / (2.0 * alpha * alpha)) / (alpha * std::sqrt(constants::PI));
                //energy /= 1 - math::erfc_x(R * std::sqrt(2.0) / (2.0 * alpha)) - R * std::sqrt(2.0) * std::exp(-R*R / (2.0 * alpha * alpha)) / (alpha * std::sqrt(constants::PI));
                //energy = alpha * std::sqrt(constants::PI) * (math::erf_x(std::sqrt(2.0) * R / (2.0 * alpha)) - math::erf_x(std::sqrt(2.0) * dist / (2.0 * alpha))) * std::exp(R*R / (2.0 * alpha*alpha)) + std::sqrt(2.0)*(dist - R);
                //energy /= dist * (math::erf_x(std::sqrt(2.0) * R / (2.0 * alpha)) * std::exp(R*R / (2.0 * alpha * alpha)) * alpha * std::sqrt(constants::PI) - R*std::sqrt(2.0));
                energy = math::erfc_x(eta * q) - math::erfc_x(eta) - (1.0 - q) * 2.0 * eta / std::sqrt(constants::PI) * std::exp(-eta * eta);
                energy /= 1.0 - math::erfc_x(eta) - 2.0 * eta / std::sqrt(constants::PI) * std::exp(-eta * eta);

                //printf("Real %.15lf\n", energy);
                return energy * q1 * q2 / dist;
            }
        }

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };









    class LongTruncated{
        private:
        std::vector<double> resFac, kNorm;
        std::vector< Eigen::Vector3d > kVec;
        std::vector< std::complex<double> > rkVec;
        double volume, selfTerm = 0.0, xb, yb, zb;

        public:

        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
        }


        void initialize(Particles &particles){
            double k2 = 0;


            printf("Setting up truncated ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            //get k-vectors
            double factor = 1;
            Eigen::Vector3d vec;
            //printf("Calculating k-vectors");
            for(int kx = -kM[0]; kx <= kM[0]; kx++){
                for(int ky = -kM[1]; ky <= kM[1]; ky++){
                    for(int kz = -kM[2]; kz <= kM[2]; kz++){
                        //if(kx^2 + ky^2+ kz^2 > Kmax^2)
                        //continue;

                        /*factor = 1.0;
                        if(kx > 0){
                            factor *= 2.0;
                        }*/

                        vec[0] = (2.0 * constants::PI * kx / this->xb);
                        vec[1] = (2.0 * constants::PI * ky / this->yb);
                        vec[2] = (2.0 * constants::PI * kz / this->zb);
                        k2 = math::dot(vec, vec);

                        if(fabs(k2) > 1e-12) {
                            if(spherical){
                                if(kx * kx + ky * ky + kz * kz < kMax * kMax){
                                    this->kVec.push_back(vec);
                                    this->resFac.push_back(factor * std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                                }
                            }
                            else{
                                this->kVec.push_back(vec);
                                this->resFac.push_back(factor * std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                            }
                        }
                    }
                }
            }

            printf("\tFound: %lu k-vectors\n", kVec.size());
            printf("\tAlpha is set to: %lf\n", alpha);
            //Calculate norms
            for(unsigned int i = 0; i < kVec.size(); i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }

            std::complex<double> rho;
            std::complex<double> rk;
            std::complex<double> charge;

            for(unsigned int k = 0; k < kVec.size(); k++){
                rho = 0;
                for(unsigned int i = 0; i < particles.tot; i++){
                    rk.imag(std::sin(math::dot(particles[i]->pos, kVec[k])));
                    rk.real(std::cos(math::dot(particles[i]->pos, kVec[k])));
                    charge = particles[i]->q;
                    rk = rk * charge;
                    rho += rk;
                }
                this->rkVec.push_back(rho);
            }

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
            //this->selfTerm *= std::sqrt(2.0) * (1.0 -  std::exp(-R*R / (2.0 * alpha * alpha)));
            this->selfTerm *= 1.0 / (std::sqrt(2.0) * alpha) / sqrt(constants::PI) * (1.0 - std::exp(-eta * eta));
            this->selfTerm /= 1.0 - math::erfc_x(eta) - 2.0 * eta / std::sqrt(constants::PI) * std::exp(-eta * eta);
            //this->selfTerm /= (1.0 - math::erfc_x(R / (std::sqrt(2.0) * alpha)) - std::sqrt(2.0) * R * std::exp(-R*R / (2.0 * alpha * alpha)) / (std::sqrt(constants::PI) * alpha)) * (std::sqrt(constants::PI) * alpha) * 2.0;
            //this->selfTerm *= alpha / sqrt(constants::PI);
            printf("\tEwald initialization Complete\n");
        }


        inline std::complex<double> Ak(unsigned int i){
            //std::complex<double> energy1;
            //std::complex<double> energy2;
            //std::complex<double> energy3;
            std::complex<double> energy4;

            std::complex<double> c1;
            std::complex<double> c2;
            std::complex<double> c3;

            double eta = R / (std::sqrt(2.0) * alpha);


            std::complex<double> zf;
            std::complex<double> zcf;
            zf.real(-kNorm[i] * R / (2.0 * eta));
            zf.imag(eta);
            zcf.real(kNorm[i] * R / (2.0 * eta));
            zcf.imag(eta);

            c1.real(R / (std::sqrt(2.0) * alpha));
            c1.imag(kNorm[i] * alpha / std::sqrt(2.0));

            c2.real(R / (std::sqrt(2.0)* alpha));
            c2.imag(-kNorm[i] * alpha / std::sqrt(2.0));

            std::complex<double> e1;
            std::complex<double> e2;
            e1.real(std::cos(kNorm[i] * R));
            e1.imag(std::sin(kNorm[i] * R));
            e2.real(std::cos(kNorm[i] * R));
            e2.imag(-std::sin(kNorm[i] * R));

            //std::cout << Faddeeva::w(zcf) << " " << Faddeeva::w(zf) << " " << (Faddeeva::w(zcf) * e1  + Faddeeva::w(zf) * e2) / 2.0 << std::endl;
            //energy1 = (Faddeeva::erf(c1) + Faddeeva::erf(c2)) * std::exp(-kNorm[i] * kNorm[i] * alpha * alpha / 2.0) / 2.0 - std::sqrt(2.0) * std::sin(R * kNorm[i]) * std::exp(-R*R / (2.0 * alpha * alpha)) / (std::sqrt(constants::PI) * alpha * kNorm[i]);
            //energy2 = (1.0 - ( (Faddeeva::w(zcf) * e1  + Faddeeva::w(zf) * e2) / 2.0 + std::sin(R * kNorm[i]) / (R * kNorm[i]) * 2.0 * eta / std::sqrt(constants::PI)) * std::exp(kNorm[i]*kNorm[i] * R*R / (4.0 * eta*eta) - eta*eta)) * std::exp(-kNorm[i]*kNorm[i] * R*R / (4.0 * eta * eta));
            //energy3 = std::exp(-kNorm[i]*kNorm[i] * R*R / (4.0 * eta * eta)) - ( (Faddeeva::w(zcf) * e1  + Faddeeva::w(zf) * e2) / 2.0 + std::sin(R * kNorm[i]) / (R * kNorm[i]) * 2.0 * eta / std::sqrt(constants::PI)) * std::exp(-eta*eta);
            c3 = Faddeeva::w(zf) * e2;
            energy4 = std::exp(-kNorm[i]*kNorm[i] * R*R / (4.0 * eta * eta)) - ( c3.real() + std::sin(R * kNorm[i]) / (R * kNorm[i]) * 2.0 * eta / std::sqrt(constants::PI)) * std::exp(-eta*eta);

            double den = 1.0 - math::erfc_x(R / (std::sqrt(2.0) * alpha)) - R * std::sqrt(2.0) * std::exp(-R*R / (2.0 * alpha * alpha)) / (std::sqrt(constants::PI) * alpha);
            //energy1 /= den;
            //energy2 /= den;
            //energy3 /= den;    
            energy4 /= den;       

            //if((std::fabs(energy1.real() - energy2.real()) > 1E-14) || (std::fabs(energy1.real() - energy3.real()) > 1E-14 )){
            //    printf("Energy difference\n");
                //exit(0);
            //}

            //printf("Ak1 real: %.15lf, complex1: %.15lf\n", energy1.real(), energy1.imag());
            //printf("Ak2 real: %.15lf, complex2: %.15lf\n", energy2.real(), energy2.imag());
            //printf("Ak3 real: %.15lf, complex3: %.15lf\n", energy3.real(), energy3.imag());
            //printf("Ak4 real: %.15lf, complex4: %.15lf\n", energy4.real(), energy4.imag());
            return energy4;
        }



        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            std::complex<double> rk_new;
            std::complex<double> rk_old;

            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto o : _old){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = o->pos.dot(this->kVec[k]);//math::dot(o->pos, this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * o->q;
                    }
                }
            }
            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto n : _new){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = n->pos.dot(this->kVec[k]);//math::dot(n->pos, this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * n->q;
                    }
                }
            }
        }


        inline double operator()(){
            double energy = 0.0;

            //#pragma omp parallel for reduction(+:energy)
            std::complex<double> _Ak;
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    _Ak = Ak(k);
                    energy += std::norm(this->rkVec[k]) * _Ak.real() * 1.0 / (kNorm[k] * kNorm[k]);//this->resFac[k];
                    if(std::fabs(_Ak.imag()) > 1E-12){
                        printf("Imaginary is too large! \n");
                        printf("%lf\n", std::fabs(_Ak.imag()));
                        exit(0);
                    }
            }
            //printf("Reciprocal term: %.15lf selfterm: %.15lf\n", energy * 2.0 * constants::PI / this->volume, this->selfTerm);
            return energy * 2.0 * constants::PI / this->volume - this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };

















    //In GC ewald should only return reciprocal part in previous state
    class Long{
        private:
        std::vector<double> resFac, kNorm;
        std::vector< Eigen::Vector3d > kVec;
        std::vector< std::complex<double> > rkVec;
        double volume, selfTerm, xb, yb, zb;

        public:

        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
            //printf("Setting volume in energy to: %lf\n", this->volume);
        }

        void set_kvectors(){
            this->kVec.clear();
            this->resFac.clear();
            this->kNorm.clear();

            double k2 = 0;

            //printf("Setting up k-vectors\n");
            //printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            //get k-vectors
            double factor = 1;
            Eigen::Vector3d vec;
            //printf("Calculating k-vectors");
            //printf("%lf %lf %lf\n", this->xb, this->yb, this->zb);
            for(int kx = -kM[0]; kx <= kM[0]; kx++){
                for(int ky = -kM[1]; ky <= kM[1]; ky++){
                    for(int kz = -kM[2]; kz <= kM[2]; kz++){

                        factor = 1.0;
                        /*if(kx > 0){
                            factor *= 2.0;
                        }*/

                        vec[0] = (2.0 * constants::PI * kx / this->xb);
                        vec[1] = (2.0 * constants::PI * ky / this->yb);
                        vec[2] = (2.0 * constants::PI * kz / this->zb);
                        k2 = math::dot(vec, vec);

                        if(fabs(k2) > 1e-12) {
                            if(spherical){
                                if(kx * kx + ky * ky + kz * kz < kMax * kMax){
                                    this->kVec.push_back(vec);
                                    this->resFac.push_back(factor * std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                                }
                            }

                            else{
                                this->kVec.push_back(vec);
                                this->resFac.push_back(factor * std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                            }
                        }
                    }
                }
            }

            //printf("\t Created %lu reciprocal lattice vectors\n", this->kVec.size());

            //Calculate norms
            for(unsigned int i = 0; i < kVec.size(); i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }
        }

        void set_self(Particles &particles){
            this->selfTerm = 0.0;

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
            this->selfTerm *= alpha / sqrt(constants::PI);
            //printf("Self term is: %lf\n", this->selfTerm);
        }


        void initialize(Particles &particles){
            set_kvectors();
            set_self(particles);
            this->rkVec.clear();

            std::complex<double> rho;
            std::complex<double> rk;
            std::complex<double> charge;

            //#pragma omp parallel for private(rk_new, rk_old)
            for(unsigned int k = 0; k < kVec.size(); k++){
                rho = 0;
                for(unsigned int i = 0; i < particles.tot; i++){
                    rk.imag(std::sin(math::dot(particles[i]->pos, kVec[k])));
                    rk.real(std::cos(math::dot(particles[i]->pos, kVec[k])));
                    charge = particles[i]->q;
                    rk = rk * charge;
                    rho += rk;
                }
                this->rkVec.push_back(rho);
            }
            //printf("\tEwald initialization Complete\n");
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            std::complex<double> rk_new;
            std::complex<double> rk_old;

            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto o : _old){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[2] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = o->pos.dot(this->kVec[k]);//math::dot(o->pos, this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * o->q;
                    }
                }
            }
            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto n : _new){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[2] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = n->pos.dot(this->kVec[k]);//math::dot(n->pos, this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * n->q;
                    }
                }
            }
        }


        inline double operator()(){
            double energy = 0.0;

            #pragma omp parallel for reduction(+:energy)
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    energy += std::norm(this->rkVec[k]) * this->resFac[k];
            }
            //printf("Reciprocal term: %.15lf selfterm: %.15lf\n", energy * 2.0 * constants::PI / (this->volume), this->selfTerm);
            return energy * 2.0 * constants::PI / (this->volume) - this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d rec = Eigen::Vector3d::Zero();
            double kSq = 0.0;

            //Reciprocal force
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                kSq = this->kNorm[k] * this->kNorm[k];
                rec += this->kVec[k] / kSq * std::exp(-kSq / (4.0 * alpha)) * std::sin(this->kVec[k].dot(disp));
            }
            rec *= 4.0 * constants::PI / this->volume;
            /*printf("\nReciprocal force: \n");
            std::cout << q1*q2*rec << std::endl;
            printf("\n");*/
            return q1 * q2 * rec;
        }
    };

















    class Long2{
        private:
        std::vector<double> resFac, kNorm;
        std::vector< Eigen::Vector3d > kVec;
        std::vector< std::complex<double> > rkVec;
        double volume, selfTerm, xb, yb, zb;

        public:

        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
            //printf("Setting volume in energy to: %lf\n", this->volume);
        }

        void set_kvectors(){
            this->kVec.clear();
            this->resFac.clear();
            this->kNorm.clear();

            double k2 = 0;

            //printf("Setting up k-vectors\n");
            //printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            //get k-vectors
            Eigen::Vector3d vec;
            //printf("Calculating k-vectors");
            //printf("%lf %lf %lf\n", this->xb, this->yb, this->zb);
            for(int kx = -kM[0]; kx <= kM[0]; kx++){
                for(int ky = -kM[1]; ky <= kM[1]; ky++){
                    vec[0] = (2.0 * constants::PI * kx / this->xb);
                    vec[1] = (2.0 * constants::PI * ky / this->yb);
                    vec[2] = 0.0;
                    k2 = math::dot(vec, vec);

                    if(fabs(k2) > 1e-10) {
                        this->kVec.push_back(vec);
                        this->resFac.push_back(std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                    }
                }
            }
            unsigned int st = this->kVec.size();
            for(int kz = -kM[2]; kz <= kM[2]; kz++){

                vec[0] = 0.0;
                vec[1] = 0.0;
                vec[2] = (2.0 * constants::PI * kz / this->zb);
                k2 = math::dot(vec, vec);

                if(fabs(k2) > 1e-10) {
                    this->kVec.push_back(vec);
                    this->resFac.push_back(std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                }
            }
            for(unsigned int i = 0; i < st; i++){
                for(int kz = -kM[2]; kz <= kM[2]; kz++){
                    if(kz == 0) continue;
                    vec[0] = this->kVec[i][0];
                    vec[1] = this->kVec[i][1];
                    vec[2] = (2.0 * constants::PI * kz / this->zb);
                    k2 = math::dot(vec, vec);

                    if(fabs(k2) > 1e-12) {
                        this->kVec.push_back(vec);
                        this->resFac.push_back(std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                    }
                }
            }
            printf("\t Created %lu reciprocal lattice vectors\n", this->kVec.size());
            //Calculate norms
            for(unsigned int i = 0; i < kVec.size(); i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }
        }

        void set_self(Particles &particles){
            this->selfTerm = 0.0;

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
            this->selfTerm *= alpha / sqrt(constants::PI);
            printf("Self term is: %lf\n", this->selfTerm);
        }


        void initialize(Particles &particles){
            set_kvectors();
            set_self(particles);
            this->rkVec.clear();

            std::complex<double> rho;
            std::complex<double> rk;
            std::complex<double> charge;

            //#pragma omp parallel for private(rk_new, rk_old)
            for(unsigned int k = 0; k < kVec.size(); k++){
                rho = 0;
                for(unsigned int i = 0; i < particles.tot; i++){
                    rk.imag(std::sin(math::dot(particles[i]->pos, kVec[k])));
                    rk.real(std::cos(math::dot(particles[i]->pos, kVec[k])));
                    charge = particles[i]->q;
                    rk = rk * charge;
                    rho += rk;
                }
                this->rkVec.push_back(rho);
            }
            //printf("\tEwald initialization Complete\n");
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            std::complex<double> rk_new;
            std::complex<double> rk_old;

            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto o : _old){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[2] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = o->pos.dot(this->kVec[k]);//math::dot(o->pos, this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * o->q;
                    }
                }
            }
            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto n : _new){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[2] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = n->pos.dot(this->kVec[k]);//math::dot(n->pos, this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * n->q;
                    }
                }
            }
        }


        inline double operator()(){
            double energy = 0.0;

            #pragma omp parallel for reduction(+:energy)
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    energy += std::norm(this->rkVec[k]) * this->resFac[k];
            }
            //printf("Reciprocal term: %.15lf selfterm: %.15lf\n", energy * 2.0 * constants::PI / (this->volume), this->selfTerm);
            return energy * 2.0 * constants::PI / (this->volume) - this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };





















    class SlabCorr{
        private:
        double dipoleMoment;
        double fac;

        public:

        void set_box(double x, double y, double z){
            this->fac = 2.0 * constants::PI / (x * y * z);
        }

        void initialize(Particles &particles){
            dipoleMoment = 0.0;
            for(unsigned int i = 0; i < particles.tot; i++){
                this->dipoleMoment += particles[i]->q * particles[i]->pos[2];
            }
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            for(auto o : _old){
                this->dipoleMoment -= o->q * o->pos[2];
            }

            for(auto n : _new){
                this->dipoleMoment += n->q * n->pos[2];
            }
        }

        inline double operator()(){
            return this->fac * this->dipoleMoment * this->dipoleMoment;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };



    class ChargeCorr{
        private:
        double volume2;
        double totCharge;

        public:

        void set_box(double x, double y, double z){
            this->volume2 = 2.0 * x * y * z;
        }

        void initialize(Particles &particles){
            this->totCharge = 0.0;
            for(unsigned int i = 0; i < particles.tot; i++){
                this->totCharge += particles[i]->q;
            }
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            for(auto o : _old){
                this->totCharge -= o->q;
            }

            for(auto n : _new){
                this->totCharge += n->q;
            }
        }

        inline double operator()(){
            return - constants::PI * this->totCharge * this->totCharge / (volume2 * alpha * alpha);
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };




    class SlabCorr2{
        private:
        Eigen::Vector3d dipoleMoment;
        double fac;

        public:

        void set_box(double x, double y, double z){
            this->fac = 2.0 * constants::PI / (3.0 * x * y * z);
        }

        void initialize(Particles &particles){
            dipoleMoment << 0.0, 0.0, 0.0;
            for(unsigned int i = 0; i < particles.tot; i++){
                this->dipoleMoment += particles[i]->q * particles[i]->pos;
            }
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            for(auto o : _old){
                this->dipoleMoment -= o->q * o->pos;
            }

            for(auto n : _new){
                this->dipoleMoment += n->q * n->pos;
            }
        }

        inline double operator()(){
            return this->fac * this->dipoleMoment.dot(this->dipoleMoment);
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };

    class SlabCorr3{
        private:
        double dipoleMoment;
        double totQ;
        double qs;
        double fac;
        double xb;
        double yb;
        double zb;
        public:

        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->fac = 2.0 * constants::PI / (x * y * z);
        }

        void initialize(Particles &particles){
            this->dipoleMoment = 0.0;
            this->totQ = 0.0;
            this->qs = 0.0;

            for(unsigned int i = 0; i < particles.tot; i++){
                this->dipoleMoment += particles[i]->q * particles[i]->pos[2];
                this->totQ += particles[i]->q;
                this->qs += particles[i]->q * particles[i]->pos[2] * particles[i]->pos[2];
            }
            printf("\tdip: %lf\n", this->dipoleMoment);
            printf("\ttotQ: %lf\n", this->totQ);
            printf("\tqs: %lf\n", this->qs);
            printf("\tb corr: %lf\n", constants::PI * this->totQ*this->totQ / (2.0 * alpha*alpha * this->xb*this->yb*this->zb));
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            for(auto o : _old){
                this->dipoleMoment -= o->q * o->pos[2];
                this->qs -= o->q * o->pos[2] * o->pos[2];
                this->totQ -= o->q;
            }

            for(auto n : _new){
                this->dipoleMoment += n->q * n->pos[2];
                this->qs += n->q * n->pos[2] * n->pos[2];
                this->totQ += n->q;
            }
        }

        inline double operator()(){
            return this->fac * (this->dipoleMoment * this->dipoleMoment - this->totQ*this->qs - this->totQ*this->totQ * this->zb*this->zb / 12.0);
        }

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };









    class LongHWIPBC{
        private:
        std::vector<double> resFac, kNorm;
        std::vector< Eigen::Vector3d > kVec;
        std::vector< std::complex<double> > rkVec;
        double volume, selfTerm = 0.0, xb, yb, zb;

        public:
        
        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
        }


        void initialize(Particles &particles){
            double k2 = 0;
            //int zMax = (int) (this->zb / this->xb * kMax);
            printf("Setting up ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            //get k-vectors
            double factor = 1;
            Eigen::Vector3d vec;
            //printf("Calculating k-vectors");
            for(int kx = 0; kx <= kM[0]; kx++){
                for(int ky = 0; ky <= kM[1]; ky++){
                    for(int kz = 0; kz <= kM[2]; kz++){

                        factor = 1.0;
                        if(kx > 0){
                            factor *= 2.0;
                        }
                        if(ky > 0){
                            factor *= 2.0;
                        }
                        if(kz > 0){
                            factor *= 2.0;
                        }

                        vec[0] = (2.0 * constants::PI * kx / this->xb);
                        vec[1] = (2.0 * constants::PI * ky / this->yb);
                        vec[2] = (2.0 * constants::PI * kz / this->zb);
                        k2 = math::dot(vec, vec);

                        if(fabs(k2) > 1e-8){// && fabs(k2) < kMax) {
                            this->kVec.push_back(vec);
                            this->resFac.push_back(factor * std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                        }
                    }
                }
            }

            printf("\tFound: %lu k-vectors\n", kVec.size());
            printf("\tAlpha is set to: %lf\n", alpha);
            //Calculate norms
            for(unsigned int i = 0; i < kVec.size(); i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }

            std::complex<double> rho;
            std::complex<double> rk;
            std::complex<double> charge;
            Eigen::Vector3d temp;

            for(unsigned int k = 0; k < kVec.size(); k++){
                rho = 0.0;
                for(unsigned int i = 0; i < particles.tot; i++){
                    
                    double cosXY = std::cos(particles[i]->pos[0] * kVec[k][0]) * std::cos(particles[i]->pos[1] * kVec[k][1]);
                    rk.imag(-cosXY * std::sin(particles[i]->pos[2] * kVec[k][2]));
                    rk.real(cosXY * std::cos(particles[i]->pos[2] * kVec[k][2]));
                    charge = particles[i]->q;
                    rk *= charge;
                    rho += rk;

                    //Mirror images
                    temp = particles[i]->pos;
                    temp[2] = math::sgn(temp[2]) * this->zb / 2.0 - temp[2]; 
                    cosXY = std::cos(temp[0] * kVec[k][0]) * std::cos(temp[1] * kVec[k][1]);
                    rk.imag(-cosXY * std::sin(temp[2] * kVec[k][2]));
                    rk.real(cosXY * std::cos(temp[2] * kVec[k][2]));
                    charge = -particles[i]->q;
                    rk *= charge;
                    rho += rk;
                }
                this->rkVec.push_back(rho);
            }

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }

            this->selfTerm *= alpha / sqrt(constants::PI); //   *2.0 due to images
            printf("\tEwald initialization Complete\n");
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            std::complex<double> rk_new;
            std::complex<double> rk_old;
            std::complex<double> charge;
            Eigen::Vector3d temp;

            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q * alpha / sqrt(constants::PI);
                }
            }
            else{
                for(auto o : _old){
                    //#pragma omp parallel for private(rk_new, rk_old)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double cosXY = std::cos(o->pos[0] * kVec[k][0]) * std::cos(o->pos[1] * kVec[k][1]);
                        rk_old.imag(-cosXY * std::sin(o->pos[2] * kVec[k][2]));
                        rk_old.real(cosXY * std::cos(o->pos[2] * kVec[k][2]));
                        charge = o->q;
                        this->rkVec[k] -= rk_old * charge;

                        // Remove image
                        temp = o->pos;
                        temp[2] = math::sgn(temp[2]) * this->zb / 2.0 - temp[2]; 
                        cosXY = std::cos(temp[0] * kVec[k][0]) * std::cos(temp[1] * kVec[k][1]);
                        rk_old.imag(-cosXY * std::sin(temp[2] * kVec[k][2]));
                        rk_old.real(cosXY * std::cos(temp[2] * kVec[k][2]));
                        charge = -o->q;
                        this->rkVec[k] -= rk_old * charge;
                    }
                }
            }
            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q * alpha / sqrt(constants::PI);
                }
            }
            else{
                for(auto n : _new){
                    //#pragma omp parallel for private(rk_new, rk_old)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double cosXY = std::cos(n->pos[0] * kVec[k][0]) * std::cos(n->pos[1] * kVec[k][1]);
                        rk_new.imag(-cosXY * std::sin(n->pos[2] * kVec[k][2]));
                        rk_new.real(cosXY * std::cos(n->pos[2] * kVec[k][2]));
                        charge = n->q;
                        this->rkVec[k] += rk_new * charge;

                        // Add image
                        temp = n->pos;
                        temp[2] = math::sgn(temp[2]) * this->zb / 2.0 - temp[2]; 
                        cosXY = std::cos(temp[0] * kVec[k][0]) * std::cos(temp[1] * kVec[k][1]);
                        rk_new.imag(-cosXY * std::sin(temp[2] * kVec[k][2]));
                        rk_new.real(cosXY * std::cos(temp[2] * kVec[k][2]));
                        charge = -n->q;
                        this->rkVec[k] += rk_new * charge;
                    }
                }
            }
        }


        inline double operator()(){
            double energy = 0.0;

            //#pragma omp parallel for reduction(+:energy)
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    energy += std::norm(this->rkVec[k]) * this->resFac[k];
            }
            
            return energy * constants::PI / (this->volume) - this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };



















    class LongHW{
        private:
        std::vector<double> resFac, kNorm;
        std::vector< Eigen::Vector3d > kVec;
        std::vector< std::complex<double> > rkVec;
        double volume, selfTerm, xb, yb, zb;

        public:
 
        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
        }



        void initialize(Particles &particles){
            double k2 = 0;

            printf("Setting up ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            this->kVec.clear();
            this->resFac.clear();
            this->kNorm.clear();
            this->rkVec.clear();
            if(!this->rkVec.empty()){
                printf("rkVec is not empty!\n");
                exit(0);
            }
            this->selfTerm = 0.0;
            
            //get k-vectors
            double factor = 1;
            //std::vector<double> vec(3);
            Eigen::Vector3d vec;
            //printf("Calculating k-vectors");
            for(int kx = 0; kx <= kM[0]; kx++){
                for(int ky = -kM[1]; ky <= kM[1]; ky++){
                    for(int kz = -kM[2]; kz <= kM[2]; kz++){

                        factor = 1.0;
                        if(kx > 0){
                            factor *= 2.0;
                        }

                        vec[0] = (2.0 * constants::PI * kx / this->xb);
                        vec[1] = (2.0 * constants::PI * ky / this->yb);
                        vec[2] = (2.0 * constants::PI * kz / this->zb);
                        k2 = math::dot(vec, vec);

                        if(fabs(k2) > 1e-8){// && fabs(k2) < kMax) {
                            this->kVec.push_back(vec);
                            this->resFac.push_back(factor * std::exp(-k2 / (4.0 * alpha * alpha)) / k2);
                        }
                    }
                }
            }

            printf("\tFound: %lu k-vectors\n", kVec.size());
            printf("\tAlpha is set to: %lf\n", alpha);
            //Calculate norms
            for(unsigned int i = 0; i < kVec.size(); i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }

            std::complex<double> rho;
            std::complex<double> rk;
            std::complex<double> charge;
            Eigen::Vector3d temp;

            for(unsigned int k = 0; k < kVec.size(); k++){
                rho = 0.0;
                for(unsigned int i = 0; i < particles.tot; i++){
                    rk.imag(std::sin(math::dot(particles[i]->pos, kVec[k])));
                    rk.real(std::cos(math::dot(particles[i]->pos, kVec[k])));
                    charge = particles[i]->q;
                    rk *= charge;
                    rho += rk;

                    //Mirror images
                    temp = particles[i]->pos;
                    temp[2] = math::sgn(temp[2]) * this->zb / 2.0 - temp[2]; 
                    rk.imag(std::sin(math::dot(temp, kVec[k])));
                    rk.real(std::cos(math::dot(temp, kVec[k])));
                    charge = -particles[i]->q;
                    rk *= charge;
                    rho += rk;
                }
                this->rkVec.push_back(rho);
            }

            

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }

            this->selfTerm *= alpha / std::sqrt(constants::PI); //   *2.0 due to images
            printf("\tSelfterm is: %lf\n", this->selfTerm);
            printf("\tEwald initialization Complete\n");
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            std::complex<double> rk_new;
            std::complex<double> rk_old;
            Eigen::Vector3d temp;

            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto o : _old){
                    temp = o->pos;
                    temp[2] = math::sgn(temp[2]) * this->zb / 2.0 - temp[2]; 

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 6)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = o->pos.dot(this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * o->q;

                        // Remove image
                        dot = temp.dot(this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * (-o->q);
                    }
                }
            }
            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto n : _new){
                    temp = n->pos;
                    temp[2] = math::sgn(temp[2]) * this->zb / 2.0 - temp[2]; 

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 6)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = n->pos.dot(this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * n->q;

                        // Add image
                        dot = temp.dot(this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * (-n->q);
                    }
                }
            }
        }


        inline double operator()(){
            double energy = 0.0;

            #pragma omp parallel for reduction(+:energy) if(kM[0] > 8)
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    energy += std::norm(this->rkVec[k]) * this->resFac[k];
            }
            return energy * constants::PI / (this->volume) - this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d rec = Eigen::Vector3d::Zero();
            double kSq = 0.0;

            //Reciprocal force
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                kSq = this->kNorm[k] * this->kNorm[k];
                rec += this->kVec[k] / (kSq) * std::exp(-kSq / (4.0 * alpha)) * 
                        std::sin(this->kVec[k].dot(disp));
            }
            rec *= 4.0 * constants::PI / this->volume;
            return q1 * q2 * rec;
        }
    };













    class LongEllipsoidal{
        private:
        std::vector<double> resFac, kNorm;
        std::vector< Eigen::Vector3d > kVec;
        std::vector< std::complex<double> > rkVec;
        double volume, selfTerm = 0.0, xb, yb, zb;

        public:

        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
        }




        void initialize(Particles &particles){
            double k2 = 0;
            double gz = 2.00000000000000;//this->xb / this->zb;
            //double cR = alpha * 5.0;
            //double cL = 5.0 / 10.0;
            //double c0 = cR / (cL * 2.0 * constants::PI);
            double c0 = (alpha * this->xb) / (2.0 * constants::PI);

            printf("Setting up ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            //get k-vectors
            double factor = 1;
            Eigen::Vector3d vec;
            Eigen::Vector3d vec2;
            //printf("Calculating k-vectors");
            for(int kx = -kM[0]; kx <= kM[0]; kx++){
                for(int ky = -kM[1]; ky <= kM[1]; ky++){
                    for(int kz = -kM[2]; kz <= kM[2]; kz++){

                        factor = 1.0;
                        //if(kx > 0){
                        //    factor *= 2.0;
                        //}

                        vec[0] = (2.0 * constants::PI * kx / this->xb);
                        vec[1] = (2.0 * constants::PI * ky / this->yb);
                        vec[2] = (2.0 * constants::PI * kz / this->zb);
                        k2 = math::dot(vec, vec);

                        vec2[0] = (kx);
                        vec2[1] = (ky);
                        vec2[2] = (kz);
                        double l2 = math::dot(vec2, vec2);

                        if(fabs(k2) > 1e-8){// && fabs(k2) < kMax) {
                            this->kVec.push_back(vec);
                            this->resFac.push_back(factor * std::exp(-l2 / (4.0 * c0 * c0)) / k2);
                        }
                    }
                }
            }

            printf("\tFound: %lu k-vectors\n", kVec.size());
            printf("\tAlpha is set to: %lf\n", alpha);
            //Calculate norms
            for(unsigned int i = 0; i < kVec.size(); i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }

            std::complex<double> rho;
            std::complex<double> rk;
            std::complex<double> charge;

            for(unsigned int k = 0; k < kVec.size(); k++){
                rho = 0;
                for(unsigned int i = 0; i < particles.tot; i++){
                    rk.imag(std::sin(math::dot(particles[i]->pos, kVec[k])));
                    rk.real(std::cos(math::dot(particles[i]->pos, kVec[k])));
                    charge = particles[i]->q;
                    rk = rk * charge;
                    rho += rk;
                }
                this->rkVec.push_back(rho);
            }

            for(unsigned int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
            //this->selfTerm *= alpha / sqrt(constants::PI);

            this->selfTerm *= (alpha * gz * (constants::PI - 2.0 * std::atan(1.0 / sqrt(gz * gz - 1.0)))) / (2.0 * sqrt(constants::PI) * sqrt(gz * gz - 1.0));
            printf("\tEwald initialization Complete\n");
        }

        inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
            std::complex<double> rk_new;
            std::complex<double> rk_old;

            if(_old.empty()){
                for(auto n : _new){
                    this->selfTerm += n->q * n->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto o : _old){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = o->pos.dot(this->kVec[k]);//math::dot(o->pos, this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * o->q;
                    }
                }
            }
            if(_new.empty()){
                for(auto o : _old){
                    this->selfTerm -= o->q * o->q * alpha / std::sqrt(constants::PI);
                }
            }
            else{
                for(auto n : _new){

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = n->pos.dot(this->kVec[k]);//math::dot(n->pos, this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * n->q;
                    }
                }
            }
        }


        inline double operator()(){
            double energy = 0.0;

            //#pragma omp parallel for reduction(+:energy)
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    energy += std::norm(this->rkVec[k]) * this->resFac[k];
            }
            printf("Reciprocal term: %.15lf selfterm: %.15lf\n", energy * 2.0 * constants::PI / (this->volume), this->selfTerm);
            return energy * 2.0 * constants::PI / (this->volume) - this->selfTerm;
        } 

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };

























    class Ewald2D{
        private:
        std::vector< Eigen::Vector3d > kVec;
        std::vector<double> kNorm;
        int kNum;
        double xb, yb, zb, volume, selfTerm;

        public:
        void set_box(double x, double y, double z){
            this->xb = x;
            this->yb = y;
            this->zb = z;
            this->volume = x * y * z;
        }

        void initialize(Particles &particles){
            kNum = 0;
            int kx = 0;
            int ky = 0;
            double k2 = 0;

            //get k-vectors
            Eigen::Vector3d vec;

            for(kx = -kM[0]; kx <= kM[0]; kx++){
                for(ky = -kM[1]; ky <= kM[1]; ky++){
                    vec[0] = (2.0 * constants::PI * kx/this->xb);
                    vec[1] = (2.0 * constants::PI * ky/this->yb);
                    vec[2] = 0.0;
                    k2 = vec.dot(vec);
                    if(fabs(k2) > 1e-5){// && ky * ky + kx * kx <= 11){
                        this->kVec.push_back(vec);
                        kNum++;
                    }
                }
            }

            printf("2D: Found: %d k-vectors\n", kNum);
            //Calculate norms
            for(int i = 0; i < kNum; i++){
                this->kNorm.push_back(math::norm(kVec[i]));
            }

            this->selfTerm = 0.0;
            for(int i = 0; i < particles.tot; i++){
                this->selfTerm += particles[i]->q * particles[i]->q;
            }
            this->selfTerm *= alpha / std::sqrt(constants::PI);
        }

        double get_self(){
            return this->selfTerm;
        }

        double f(double norm, double zDist){
            double f = std::exp(norm * zDist) * math::erfc_x(alpha * zDist + norm/(2.0 * alpha)) +
                       std::exp(-norm * zDist) * math::erfc_x(-alpha * zDist + norm/(2.0 * alpha));

            return f / norm;
        }

        double gE(double q1, double q2, Eigen::Vector3d dispVec){
            double zDist = std::abs(dispVec[2]);
            double gt = zDist * math::erf_x(alpha * zDist) + std::exp( -(zDist * zDist * alpha * alpha) ) / ( alpha * std::sqrt(constants::PI) );
            return q1 * q2 * constants::PI/(this->xb * this->yb) * gt;
        }

        double get_reciprocal(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> &p2, Eigen::Vector3d dispVec){
            double energy = 0.0;
            double rk;
            double zDist = std::abs(dispVec[2]);
        
            for(int i = 0; i < kNum; i++){
                rk = std::cos(this->kVec[i].dot(dispVec));
                energy += rk * this->f(this->kNorm[i], zDist);
            }
            return energy;
        }

        double operator()(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, double dist){
            return p1->q * p2->q * math::erfc_x(alpha * dist) / dist;
        }

        double rec(std::shared_ptr<Particle>  p1, std::shared_ptr<Particle>  p2, Eigen::Vector3d dispVec){
            double zDist = std::abs(dispVec[2]);
            double reciprocal = this->get_reciprocal(p1, p2, dispVec);
            return p1->q * p2->q * constants::PI/(2.0 * this->xb * this->yb) * reciprocal;
        }
/*
        double get_energy(Particles &particles){
            double energy = 0;
            double distance = 0;
            double real = 0;
            double self = 0;
            double reciprocal = 0;
            double reci = 0;
            double dipCorr = 0;
            double done = 0;
            double gE = 0;
            int k = 0;

            for(int i = 0; i < particles.tot; i++){

                k = i + 1;
                while(k < particles.tot){
                    distance = particles.distances[i][k];
                    if(k != i){

                        real += particles[i]->q * particles[k]->q * math::erfc_x(alpha * distance) / distance;
                    }
                    k++;
                }

                for(int j = 0; j < particles.tot; j++){
                    double zDist = particles[i]->pos[2] - particles[j]->pos[2];
                    reciprocal += particles[i]->q * particles[j]->q * get_reciprocal(particles[i], particles[j]) * 1.0 / 2.0;
                    gE += particles[i]->q * particles[j]->q * (zDist * math::erf_x(alpha * zDist) +
                        exp(-(zDist * zDist * alpha * alpha)) / (alpha * sqrt(constants::PI)));
                }
                //dipCorr += dipole_correction(particles[i]);
                self += get_self_correction(particles[i]);
                printf("Rec: %lf g: %lf: reci: %lf Done: %lf\r", reciprocal, gE, reci,(double)i/particles.tot * 100);
            }
            reciprocal = (reciprocal - gE) * constants::PI/(Base::xL * Base::yL);
            dipCorr = -2.0 * constants::PI/(this->xb * this->yb * this->zb) * dipCorr * dipCorr;
            self = alpha/std::sqrt(constants::PI) * self;
            energy = (real + reciprocal - self);// + dipCorr);
            printf("Real: %lf, self: %lf, reciprocal: %lf dipCorr: %lf, energy: %lf\n", real * Base::lB, self * Base::lB, reciprocal * Base::lB, dipCorr, energy * Base::lB);
            return energy * constants::lB;
        }
        */

        inline Eigen::Vector3d force(double q1, double q2, Eigen::Vector3d disp){
            Eigen::Vector3d force;
            return force;
        }
    };
}


















/*

class Levin{

    private:
    double sumC, xb, yb, zb, volume, dipoleM;
    std::vector<double> f1, f2, f3, f4, kNorms;
    std::vector<int> kM;
    std::vector< Eigen::Vector2d> kVec;

    public:

    void set_km(std::vector<int> v){
        kM = v;
    }

    void set_box(double x, double y, double z){
        this->xb = x;
        this->yb = y;
        this->zb = z;
        this->volume = x * y * z;
    }

    void initialize(Particles &particles){
        Eigen::Vector2d vec;

        for(int x = -kM[0]; x <= kM[0]; x++){
            for(int y = -kM[1]; y <= kM[1]; y++){
                if(x != 0 || y != 0){
                    vec[0] = (double) x;//(double) x * 2.0 * PI / Base::xL; // 
                    vec[1] = (double) y;//(double) y * 2.0 * PI / Base::xL; // 
                    kVec.push_back(vec);
                    kNorms.push_back(2.0 * PI * sqrt((double)(x * x)/(Base::xL * Base::xL) + (double)(y * y)/(Base::yL * Base::yL)));

                }
            }
        }
        printf("\tFound %lu k-vectors.\n", kVec.size());

        double factor = 0;
        for(int i = 0; i < kVec.size(); i++){
            for(int j = 0; j < particles.tot; j++){
                factor = 2.0 * constants::PI / Base::xL * (kVec[i][0] * particles[j]->pos[0] + kVec[i][1] * particles[j]->pos[1]);
                double fac = kNorms[i] * (particles[j]->pos[2] + this->zb / 2.0);
                f1[i] += particles[j]->q * std::cos(factor) * std::exp(-fac);
                f2[i] += particles[j]->q * std::sin(factor) * std::exp(-fac);
                f3[i] += particles[j]->q * std::cos(factor) * std::exp( fac);
                f4[i] += particles[j]->q * std::sin(factor) * std::exp( fac);
            }  
        }

        printf("Calculated f-functions\n");
        for(int i = 0; i < kVec.size(); i++){
            eFactors[i] = exp(-2.0 * kNorms[i] * Base::zLBox);
        }

        for(int i = 0; i < particles.tot; i++){
            sumC += particles[i]->q;
        }

        for(int i = 0; i < particles.numOfParticles; i++){
            dipoleM += particles[i].q * ( particles[i].pos[2] + this->zb / 2.0 );
        }
    }

    inline double operator()(){
        double polarization = 0.0, gamma = 0.0, dipoleM = 0.0;;

        for(int i = 0; i < kVec.size(); i++){    
            polarization += -1.0 / (kNorms[i] * (1.0 - eFactors[i])) * 
                    (f1[i] * f1[i] + f2[i] * f2[i] + eFactors[i] * (f3[i] * f3[i] + f4[i] * f4[i]) - 
                    2.0 * eFactors[i] * (f3[i] * f1[i] + f2[i] * f4[i]));
        }

        gamma = -2.0 * (dipoleM * dipoleM / this->zb - sumC * dipoleM);

        return constants::PI / (Base::xL * Base::xL) * (polarization + gamma);
    }

    void update_f(Particle &_old, Particle &_new){
        for(int i = 0; i < kVec.size(); i++){
            double oldFactor = 2.0 * constants::PI / this->xb * (kVec[i][0] * _old.pos[0] + kVectors[i][1] * _old.pos[1]);
            double newFactor = 2.0 * constants::PI / this->xb * (kVec[i][0] * _new.pos[0] + kVectors[i][1] * _new.pos[1]);
            f1[i] -= _old.q * std::cos(oldFactor) * std::exp(-kNorms[i] * (_old.pos[2] + this->zb/ 2.0));
            f1[i] += _new.q * std::cos(newFactor) * std::exp(-kNorms[i] * (_new.pos[2] + this->zb / 2.0));

            f2[i] -= _old.q * std::sin(oldFactor) * std::exp(-kNorms[i] * (_old.pos[2] + this->zb / 2.0));
            f2[i] += _new.q * std::sin(newFactor) * std::exp(-kNorms[i] * (_new.pos[2] + this->zb / 2.0) );

            f3[i] -= _old.q * std::cos(oldFactor) * std::exp(kNorms[i] * (_old.pos[2] + this->zb / 2.0) );
            f3[i] += _new.q * std::cos(newFactor) * std::exp(kNorms[i] * (_new.pos[2] + this->zb / 2.0) );

            f4[i] -= _old.q * std::sin(oldFactor) * std::exp(kNorms[i] * (_old.pos[2] + this->zb / 2.0) );
            f4[i] += _new.q * std::sin(newFactor) * std::exp(kNorms[i] * (_new.pos[2] + this->zb / 2.0));

            sumC += _new.q;
            sumC -= _old.q;

            dipoleM += 
        }
    }
};

*/
