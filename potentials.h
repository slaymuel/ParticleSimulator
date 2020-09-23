#include "particle.h"
#include <vector>
#include "geometry.h"
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
        if(dist * dist > this->Rsq) return 1e30;
        return this->k * this->Rsq * (this->Rsq / (this->Rsq - dist * dist) - 1.0);
    }
};



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

            //printf("\tFound: %lu k-vectors\n", kVec.size());
            //printf("\tAlpha is set to: %lf\n", alpha);
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

            #pragma omp parallel for reduction(+:energy)
            for(unsigned int k = 0; k < this->kVec.size(); k++){
                    energy += std::norm(this->rkVec[k]) * this->resFac[k];
            }
            //printf("Reciprocal term: %.15lf selfterm: %.15lf\n", energy * 2.0 * constants::PI / (this->volume), this->selfTerm);
            return energy * 2.0 * constants::PI / (this->volume) - this->selfTerm;
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
