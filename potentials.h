#include "particle.h"
#include <vector>
#include "geometry.h"
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
    int kMax = 0.0;
    std::vector<int> kM;




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



    //In GC ewald should only return reciprocal part in previous state
    class Long{
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


            printf("Setting up ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

            //get k-vectors
            double factor = 1;
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
            this->selfTerm *= alpha / sqrt(constants::PI);
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
            int zMax = (int) (this->zb / this->xb * kMax);
            printf("Setting up ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kMax, kMax, zMax);

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

            printf("Setting up ewald\n");
            printf("\tWavevectors in x, y, z: %i, %i, %i\n", kM[0], kM[1], kM[2]);

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

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = o->pos.dot(this->kVec[k]);//math::dot(o->pos, this->kVec[k]);
                        rk_old.imag(std::sin(dot));
                        rk_old.real(std::cos(dot));

                        this->rkVec[k] -= rk_old * o->q;

                        // Remove image
                        dot = temp.dot(this->kVec[k]);//math::dot(temp, this->kVec[k]);
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

                    #pragma omp parallel for private(rk_new, rk_old) if(kM[0] > 8)
                    for(unsigned int k = 0; k < kVec.size(); k++){
                        double dot = n->pos.dot(this->kVec[k]);//math::dot(n->pos, this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * n->q;

                        // Add image
                        dot = temp.dot(this->kVec[k]);//math::dot(temp, this->kVec[k]);
                        rk_new.imag(std::sin(dot));
                        rk_new.real(std::cos(dot));

                        this->rkVec[k] += rk_new * (-n->q);
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