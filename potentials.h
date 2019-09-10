#include "particle.h"
#include <vector>
#include "geometry.h"

#pragma omp declare reduction(vec_double_plus : std::vector<std::complex<double>> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<std::complex<double>>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

class Coulomb{
    public:

    inline double operator()(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, double dist){
        return p1->q * p2->q / dist;
    }
};


class EwaldShort{
    double alpha;

    public:

    inline double operator()(std::shared_ptr<Particle>& p1, std::shared_ptr<Particle>& p2, double dist){
        this->alpha = 8.0 / 10.0;
        double energy = math::erfc_x(dist * this->alpha) / dist;
        double real = p1->q * p2->q * energy;


        return real;    //tinfoil
    }
};



//In GC ewald should only return reciprocal part in previous state
class EwaldLong{
    private:
    std::vector<double> resFac, kNorm;
    std::vector< std::vector<double> > kVec;
    std::vector< std::complex<double> > rkVec;
    double volume, alpha, selfTerm = 0.0, xb, yb, zb;

    public:

    void set_box(double x, double y, double z){
        this->xb = x;
        this->yb = y;
        this->zb = z;
        this->volume = x * y * z;
    }


    void initialize(Particles &particles){
        double k2 = 0;
        int kMax = 4;//8/this->xb;
        int zMax = (int) (this->zb / this->xb * kMax);
        this->alpha = 8.0 / 10.0;
        printf("Setting up ewald\n");
        printf("\tWavevectors in x, y, z: %i, %i, %i\n", kMax, kMax, zMax);

        //get k-vectors
        double factor = 1;
        std::vector<double> vec(3);
        //printf("Calculating k-vectors");
        for(int kx = 0; kx <= kMax; kx++){
            for(int ky = -kMax; ky <= kMax; ky++){
                for(int kz = -zMax; kz <= zMax; kz++){

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
                        this->resFac.push_back(factor * std::exp(-k2 / (4.0 * this->alpha * this->alpha)) / k2);
                    }
                }
            }
        }

        printf("\tFound: %lu k-vectors\n", kVec.size());

        //Calculate norms
        for(int i = 0; i < kVec.size(); i++){
            this->kNorm.push_back(math::norm(kVec[i]));
        }

        std::complex<double> rho;
        std::complex<double> rk;
        std::complex<double> charge;

        for(int k = 0; k < kVec.size(); k++){
            rho = 0;
            for(int i = 0; i < particles.tot; i++){
                rk.imag(std::sin(math::dot(particles[i]->pos, kVec[k])));
                rk.real(std::cos(math::dot(particles[i]->pos, kVec[k])));
                charge = particles[i]->q;
                rk = rk * charge;
                rho += rk;
            }
            this->rkVec.push_back(rho);
        }

        for(int i = 0; i < particles.tot; i++){
            this->selfTerm += particles[i]->q * particles[i]->q;
        }
        this->selfTerm *= this->alpha/sqrt(constants::PI);
        printf("\tEwald initialization Complete\n");
    }

    inline void update(std::vector< std::shared_ptr<Particle> >& _old, std::vector< std::shared_ptr<Particle> >& _new){
        std::complex<double> rk_new;
        std::complex<double> rk_old;
        std::complex<double> charge;
        //charge = 0;
        //std::cout << "rkvec before " << std::accumulate(rkVec.begin(), rkVec.end(), charge) << "\n";
        for(auto o : _old){
            #pragma omp parallel for private(rk_new, rk_old, charge)
            for(int k = 0; k < kVec.size(); k++){
                rk_old.imag(std::sin(math::dot(o->pos, this->kVec[k])));
                rk_old.real(std::cos(math::dot(o->pos, this->kVec[k])));
                charge = o->q;

                this->rkVec[k] -= rk_old * charge;
            }
        }

        for(auto n : _new){
            #pragma omp parallel for private(rk_new, rk_old, charge)
            for(int k = 0; k < kVec.size(); k++){
                //std::cout << "element before " << this->rkVec[k] << "\n";
                rk_new.imag(std::sin(math::dot(n->pos, this->kVec[k])));
                rk_new.real(std::cos(math::dot(n->pos, this->kVec[k])));
                charge = n->q;
                
                this->rkVec[k] += rk_new * charge;
                //std::cout << "element after " << this->rkVec[k] << "\n";
            }
        }
        //charge = 0;
        //std::cout << "rkvec after " << std::accumulate(this->rkVec.begin(), this->rkVec.end(), charge) << "\n";
    }


    inline double operator()(){
        double energy = 0.0;

        #pragma omp parallel for reduction(+:energy)
        for(int k = 0; k < this->kVec.size(); k++){
                energy += std::norm(this->rkVec[k]) * this->resFac[k];
        }
        
        //std::cout << "prefactors potential: " << 2.0 * constants::PI / (this->volume) << "\n";
        //std::cout << "selfTerm: " << selfTerm << "\n";
        //printf("Energy in potential: %lf\n", energy * 2.0 * constants::PI / (this->volume) );
        return energy * 2.0 * constants::PI / (this->volume) - this->selfTerm;
    } 
};