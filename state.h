#pragma once

#include <vector>
#include "particles.h"
#include <limits>
//#include <math.h>
#include "geometry.h"
#include "energy.h"
#include "potentials.h"
#include "Spline.h"

class State{
    private:

    
    SplineData spline;
    //IO io;
    
    public:

    std::shared_ptr<State> _old;
    int step = 0;
    double energy = 0.0, cummulativeEnergy = 0.0, dE = 0.0, error = 0.0;
    Particles particles;
    std::vector< unsigned int > movedParticles;    //Particles that has moved from previous state
    Geometry *geo;
    std::vector< std::shared_ptr<EnergyBase> > energyFunc;

    ~State(){
        delete geo;
    }

    void advance(){
        this->step++;
    }

    void control(){
        #ifdef DEBUG
        printf("Control (DEBUG)\n");
        #else
        printf("Control\n");
        #endif
        
        this->energy = 0.0;
        for(auto e : this->energyFunc){
            //printf("Before initialize: %.10lf\n", e->all2all(this->particles));
            if(this->step % 100 == 0){
                e->initialize(this->particles);
            }
            //printf("After initialize: %.10lf\n", e->all2all(this->particles));
            this->energy += e->all2all(this->particles);
        }

        this->error = std::fabs((this->energy - this->cummulativeEnergy) / this->energy);


        #ifdef DEBUG
        if(this->particles.tot != _old->particles.tot){
            printf("_old state has %i particles and current state has %i.\n", _old->particles.tot, this->particles.tot);
            exit(1);
        }
        unsigned int cations = 0, anions = 0;
        for(unsigned int i = 0; i < this->particles.tot; i++){
            if(std::abs(this->geo->distance(this->particles.particles[i]->com, this->particles.particles[i]->pos) - 
                                                                this->particles.particles[i]->b) > 1e-5){
                printf("|Pos - com| is not b! it is: %.6lf and should be: %.6lf\n",
                                    this->geo->distance(this->particles.particles[i]->com, this->particles.particles[i]->pos), 
                                    this->particles.particles[i]->b);
                exit(1);
            }
            if(std::abs(this->geo->distance(this->particles.particles[i]->com, this->particles.particles[i]->pos) - 
                                                                this->particles.particles[i]->qDisp.norm()) > 1e-5){
                printf("|Pos - com| is not equal to |qDisp|!\n");
                exit(1);
            }
            if(this->particles.particles[i]->b > this->particles.particles[i]->b_max){
                printf("Ooops, b is larger than b_max!\n");
                exit(1);
            }
            if(this->particles.particles[i]->pos != this->_old->particles.particles[i]->pos){
                printf("current positions is not equal to old positions for particle %i.\n", i);
                std::cout << this->particles.particles[i]->pos << std::endl;
                std::cout << "\n" << this->_old->particles.particles[i]->pos<< std::endl;
                exit(1);
            }
            if(this->particles.particles[i]->com != this->_old->particles.particles[i]->com){
                printf("current center of mass is not equal to old for particle %i.\n", i);
                std::cout << this->particles.particles[i]->com << std::endl;
                std::cout << "\n" << this->_old->particles.particles[i]->com << std::endl;
                exit(1);
            }
            if(this->particles.particles[i]->index != i){
                printf("index is wrong in current for particle %i, it has index %i.\n", i, this->particles.particles[i]->index );
                exit(1);
            }
            if(this->_old->particles.particles[i]->index != i){
                printf("index is wrong in in _old for particle %i, it has index %i.\n", i, this->_old->particles.particles[i]->index );
                exit(1);
            }

            (this->particles.particles[i]->q > 0.0) ? cations++ : anions++;
        }

        if(cations != this->particles.cTot || anions != this->particles.aTot){
            printf("Wrong number of cations or anions!\n");
            printf("Cations %i should be %i. Anions %i should be %i\n", this->particles.cTot, cations, this->particles.aTot, anions);
            exit(0);
        }

        if(this->particles.cTot + this->particles.aTot != this->particles.tot){
            printf("cTot + aTot != tot\n");
            exit(0);
        }

        if(!this->movedParticles.empty()){
            printf("Moved particles is not empty!\n");
            exit(0);
        }

        if(!this->_old->movedParticles.empty()){
            printf("_old Moved particles is not empty!\n");
            exit(0);
        }
        #endif
        if(this->error > 1e-10 || this->energy > 1e30){
            printf("\n\nEnergy drift is too large: %.12lf (all2all: %lf, cummulative: %lf)\n\n", this->error, this->energy, this->cummulativeEnergy);
            exit(1);
        } 
        /*if(this->energy != 0 && this->cummulativeEnergy != 0){

        }*/
    }

    void finalize(std::string name){
        printf("\nFinalizing simulation: %s.\n", name.c_str());
        // Set up old system
        for(std::shared_ptr<Particle> p : this->particles.particles){
            _old->particles.add(p);
        }

        //Calculate the initial energy of the system
        for(auto e : this->energyFunc){
            e->initialize(particles);
            printf("energy: %lf\n", e->all2all(this->particles));
            this->energy += e->all2all(this->particles);
        }
        this->cummulativeEnergy = this->energy;
        printf("\tEnergy of the first frame is: %.15lf\n\n", this->energy);
    }

    void reset_energy(){
        this->energy = 0.0;

        for(auto e : this->energyFunc){
            this->energy += e->all2all(this->particles);
        }

        this->cummulativeEnergy = this->energy;
    }

    void save(){
        for(auto i : this->movedParticles){
            if(this->particles.tot > this->_old->particles.tot){
                this->_old->particles.add(this->particles.particles[i]);
            }
            else if(this->particles.tot == this->_old->particles.tot){
                *(this->_old->particles.particles[i]) = *(this->particles.particles[i]);
                //For SingleSwap move
                this->_old->particles.cTot = this->particles.cTot;
                this->_old->particles.aTot = this->particles.aTot;
            }
        }
        //std::sort(this->_old->movedParticles.begin(), this->_old->movedParticles.end());
        //std::reverse(this->_old->movedParticles.begin(), this->_old->movedParticles.end());
        for(auto i : this->_old->movedParticles){
            if(this->particles.tot < this->_old->particles.tot){
                //printf("removing %i from old\n", i);
                this->_old->particles.remove(i);
            }
        }

        this->movedParticles.clear();
        this->_old->movedParticles.clear();
        this->cummulativeEnergy += this->dE;

        if(this->geo->volume != this->_old->geo->volume){
            //Update old geometry
            this->_old->geo->d = this->geo->d;
            this->_old->geo->_d = this->geo->_d;
            this->_old->geo->dh = this->geo->dh;
            this->_old->geo->_dh = this->geo->_dh;
            this->_old->geo->volume = this->geo->volume;
        }
    }


    void revert(){
        //Set moved partiles in current state equal to previous state
        if(this->dE != std::numeric_limits<double>::infinity()){
            for(auto e : this->energyFunc){
                if(this->geo->volume != this->_old->geo->volume){
                    // Is update really needed here?
                    e->update(this->_old->geo->d[0], this->_old->geo->d[1], this->_old->geo->d[2]);
                    e->initialize(this->_old->particles);
                }
                else{
                    e->update( this->particles.get_subset(this->movedParticles), this->_old->particles.get_subset(this->_old->movedParticles) );
                }
            }
        }

        //printf("p1 %.8lf %.8lf %.8lf\n", this->particles[0]->com[0], this->particles[0]->com[1], this->particles[0]->com[2]);


        if(this->particles.tot == this->_old->particles.tot){
            for(auto i : this->movedParticles){
                *(this->particles.particles[i]) = *(this->_old->particles.particles[i]);
                //For SingleSwap move
                this->particles.cTot = this->_old->particles.cTot;
                this->particles.aTot = this->_old->particles.aTot;
            }
        }
        else if(this->particles.tot > this->_old->particles.tot){   //Added particle
            std::reverse(this->movedParticles.begin(), this->movedParticles.end());
            for(auto i : this->movedParticles){
                this->particles.remove(i);
            }
        }
        else if(this->particles.tot < this->_old->particles.tot){   //Removed particle
            //std::sort(this->_old->movedParticles.begin(), this->_old->movedParticles.end());
            std::reverse(this->_old->movedParticles.begin(), this->_old->movedParticles.end());

            for(auto i : this->_old->movedParticles){
                //printf("adding back particle %i\n", i);
                this->particles.add(this->_old->particles.particles[i], i);
            }
        }
        /*for(auto i : this->movedParticles){
            if(this->particles.tot > _old->particles.tot){
                //printf("\nReject: removing particle %i from current\n\n", i);
                this->particles.remove(i);
            }
            else if(this->particles.tot == this->_old->particles.tot){
                //printf("\nAssuming normal move\n");
                *(this->particles.particles[i]) = *(this->_old->particles.particles[i]);
                //For SingleSwap move
                this->particles.cTot = this->_old->particles.cTot;
                this->particles.aTot = this->_old->particles.aTot;
            }
        }
    
        for(auto i : this->_old->movedParticles){
            if(this->particles.tot < this->_old->particles.tot){
                //printf("\nRevert: adding back particle %i to current\n", i);
                this->particles.add(this->_old->particles.particles[i], i);
            }

        }*/

        this->movedParticles.clear();
        this->_old->movedParticles.clear();

        //Revert geometry
        this->geo->d      = this->_old->geo->d;
        this->geo->_d     = this->_old->geo->_d;
        this->geo->dh     = this->_old->geo->dh;
        this->geo->_dh    = this->_old->geo->_dh;
        this->geo->volume = this->_old->geo->volume;
    }


    //Get energy different between *this and old state
    double get_energy_change(){
        double E1 = 0.0, E2 = 0.0;
        //auto start = std::chrono::steady_clock::now();
        for(auto p : this->movedParticles){
            if(!this->geo->is_inside(this->particles.particles[p]) || this->overlap(p)){
                this->dE = std::numeric_limits<double>::infinity();
                return this->dE;
            }
        }
        int counter = 0;
        for(auto e : this->energyFunc){

            /*if(this->particles.tot > this->_old->particles.tot){
                if(counter == 2){
                    continue;
                }
            }*/

            //stupid design
            e->geo = this->_old->geo;
            //std::cout << "old " << counter << " " << (*e)( this->_old->movedParticles, this->_old->particles ) << std::endl;
            E1 += (*e)( this->_old->movedParticles, this->_old->particles );

            e->geo = this->geo;
            if(this->geo->volume != this->_old->geo->volume){
                e->update(this->geo->d[0], this->geo->d[1], this->geo->d[2]);
                e->initialize(particles);
            }
            else{
                e->update( this->_old->particles.get_subset(this->_old->movedParticles), this->particles.get_subset(this->movedParticles) );
            }
            //std::cout << "new " << counter << " " << (*e)( this->movedParticles, this->particles ) << std::endl;
            E2 += (*e)( this->movedParticles, this->particles );

            counter++;
        }

        this->dE = E2 - E1;
        return this->dE;
    }


    //Called after move - set movedParticles
    void move_callback(std::vector< unsigned int > ps){  
        //this->movedParticles.insert(std::end(movedParticles), std::begin(ps), std::end(ps));

        //If a particle is removed, this->movedparticles is empty. 
        //If particle is added this->_old->particles is empty
        if(this->particles.tot >= this->_old->particles.tot){
            std::for_each(std::begin(ps), std::end(ps), [this](int i){ 
                                                    this->movedParticles.push_back(i); });
        }
        
        std::copy_if(ps.begin(), ps.end(), std::back_inserter(this->_old->movedParticles), 
                                            [this](unsigned int i){ return i < this->_old->particles.tot; });

        for(auto p : this->movedParticles){
            geo->pbc(this->particles[p]);
        }
    }


    void equilibrate(double step){
        printf("\nEquilibrating:\n");
        Eigen::Vector3d v;
        
        // Initial Check
        int i = 0, overlaps = this->get_overlaps();
        if(overlaps > 0){
            printf("\tRandomly placing particles\n");
            for(unsigned int i = 0; i < this->particles.tot; i++){
                this->particles.particles[i]->com = this->geo->random_pos(this->particles.particles[i]->rf);
                this->particles.particles[i]->pos = this->particles.particles[i]->com + this->particles.particles[i]->qDisp;
            }
        }

        printf("\tInitial overlaps: %i\n", overlaps);
        Eigen::Vector3d oldCom;
        Eigen::Vector3d oldPos;
        std::shared_ptr<Particle> p;
        double step_rand;
        //Move particles to prevent overlap
        if(overlaps > 0){
            printf("\tRemoving overlaps.\n");
            while(overlaps > 0){
                p = this->particles.random();
                oldCom = p->com;
                oldPos = p->pos;
                step_rand = Random::get_random() * step;
                p->translate(step_rand);
                //this->geo->pbc(p);
                if(!this->geo->is_inside(p) || this->overlap(p->index)){
                    p->com = oldCom;
                    p->pos = oldPos;
                }


                if(i % 50000 == 0){
                    overlaps = this->get_overlaps();
                    printf("\tOverlaps: %i, iteration: %i\r", overlaps, i);
                    fflush(stdout);
                }
                i++;
                if(i > 1E9){
                    i = 0;
                }
            }
        }
        else{
            printf("\tNo overlaps to remove!\n");
        }
        printf("\n\tEquilibration done\n\n");
    }



    bool overlap(std::size_t i){
        for(auto p : this->particles.particles){
            if(p->index == i) continue;

            if(this->geo->distance(p->com, this->particles.particles[i]->com) <= p->r + this->particles.particles[i]->r){
                return true;
            }
        }
        return false;
    }

    int get_overlaps(){
        int count = 0;
        for(auto p : this->particles.particles){
            (this->overlap(p->index)) ? count++ : 0;
        }
        return count;
    }

    void set_geometry(int type, std::vector<double> args){
        this->_old = std::make_shared<State>();

        switch (type){
            default:
                printf("Creating Cuboid box\n");
                assert(args.size() == 3);
                this->geo = new Cuboid<true, true, true>(args[0], args[1], args[2]);
                this->_old->geo = new Cuboid<true, true, true>(args[0], args[1], args[2]);
                break;

            case 1:
                this->geo = new Sphere();
                break;

            case 2:
                printf("Creating Cuboid-Image box\n");
                assert(args.size() == 3);
                this->geo = new CuboidImg<true, true, true>(args[0], args[1], args[2]);
                this->_old->geo = new CuboidImg<true, true, true>(args[0], args[1], args[2]);
                break;

            case 3:
                printf("Creating Cuboid-Image box with no PBC in z\n");
                assert(args.size() == 3);
                this->geo = new CuboidImg<true, true, false>(args[0], args[1], args[2]);
                this->_old->geo = new CuboidImg<true, true, false>(args[0], args[1], args[2]);
                break;

            case 4:
                printf("Creating Cuboid box with no PBC\n");
                assert(args.size() == 3);
                this->geo = new Cuboid<false, false, false>(args[0], args[1], args[2]);
                this->_old->geo = new Cuboid<false, false, false>(args[0], args[1], args[2]);
                break;

            case 5:
                printf("Creating Cuboid box with no PBC in z\n");
                assert(args.size() == 3);
                this->geo = new Cuboid<true, true, false>(args[0], args[1], args[2]);
                this->_old->geo = new Cuboid<true, true, false>(args[0], args[1], args[2]);
                break;
        }

    }
    


    void set_energy(int type, std::vector<double> args = std::vector<double>()){
        switch (type){
            case 1:
                printf("\nAdding Ewald potential\n");
                assert(args.size() == 7);
                this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::Short> >() );
                //this->energyFunc.push_back( std::make_shared< PairEnergyWithRep<EwaldLike::Short> >(1) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::Long> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
                EwaldLike::alpha = args[4];
                EwaldLike::kMax = args[5];
                EwaldLike::spherical = bool(args[6]);

                printf("\tSpherical cutoff: %s", EwaldLike::spherical ? "true\n" : "false\n");
                printf("\tReciprocal cutoff: %lf\n", EwaldLike::kMax);
                printf("\tk-vectors: %d %d %d\n", (int) args[1], (int) args[2], (int) args[3]);
                break;

            case 2:
                printf("\nAdding Halfwald potential\n");
                assert(args.size() == 5);
                this->energyFunc.push_back( std::make_shared< ImgEnergy<EwaldLike::Short> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongHW> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
                EwaldLike::alpha = args[4];
                printf("\tk-vectors: %d %d %d\n", (int) args[1], (int) args[2], (int) args[3]);
                break;
            
            case 3:
                printf("\nAdding HalfwaldIPBC potential\n");
                assert(args.size() == 5);
                this->energyFunc.push_back( std::make_shared< ImgEnergy<EwaldLike::Short> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongHWIPBC> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
                EwaldLike::alpha = args[4];
                break;

            case 4:
                printf("\nAdding Ellipsoidal Ewald\n");
                assert(args.size() == 5);
                this->energyFunc.push_back( std::make_shared< Ellipsoid<BSpline2D> >(spline.aKnots, spline.bKnots, spline.controlPoints) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongEllipsoidal> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
                EwaldLike::alpha = args[4];
                break;

            case 5:
                printf("\nAdding Minimum Image Halfwald\n");
                assert(args.size() == 3);
                this->energyFunc.push_back( std::make_shared< MIHalfwald<Coulomb> >(args[1], args[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                // Set box size due to reflections, needed for distance PBC
                printf("\tResetting box size in z to %lf\n", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
                this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 

                break;


            case 6:
                printf("\nAdding Truncated Ewald potential\n");
                assert(args.size() == 8);
                //this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::ShortTruncated> >() );
                this->energyFunc.push_back( std::make_shared< PairEnergyWithRep<EwaldLike::ShortTruncated> >(1) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongTruncated> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
                EwaldLike::alpha = args[4];
                //printf("Sigma: %lf\n", args[4]);
                EwaldLike::R = args[5];
                EwaldLike::kMax = args[6];
                EwaldLike::spherical = bool(args[7]);
                printf("\tSpherical cutoff: %s", EwaldLike::spherical ? "true\n" : "false\n");
                printf("\tReciprocal cutoff: %lf\n", EwaldLike::kMax);
                EwaldLike::eta = EwaldLike::R * 1.0 / (std::sqrt(2.0) * EwaldLike::alpha);
                break;

            case 7:
                printf("\nAdding Halfwald with real replicates\n");
                assert(args.size() == 7);
                this->energyFunc.push_back( std::make_shared< MIHalfwald<EwaldLike::Short> >(args[1], args[6]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                printf("\tResetting box size in z to %lf\n", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
                this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongHW> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2] * 2.0) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
                EwaldLike::alpha = args[5];

                break;

            case 8:
                printf("\nAdding harmonic well to charges\n");
                assert(args.size() == 1);
                this->energyFunc.push_back( std::make_shared< ChargeWell<Harmonic> >(args[0], 0.0) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(100.0);
                break;

            case 9:
                printf("\nAdding Sture-potential to charges\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< ChargeWell<Sture> >(args[0], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(100.0);
                break;

            case 10:
                printf("\nAdding FENE-potential to charges\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< ChargeWell<FENE> >(args[0], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(100.0);
                break;

            case 11:
                printf("\nAdding FanourgakisSP2 with image charges\n");
                assert(args.size() == 2);
                                                                                          // kMax     eps
                this->energyFunc.push_back( std::make_shared< MIHalfwald<Fanourgakis::SP2> >(args[1], 1.0) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<Fanourgakis::SP2Self> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2] * 2.0) );
                this->energyFunc.back()->set_geo(this->geo);

                Fanourgakis::R = args[0];

                printf("\tResetting box size in z to %lf\n", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
                this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 
                break;

            case 12:
                printf("\nAdding FanourgakisSP3 with image charges\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< MIHalfwald<Fanourgakis::SP3> >(args[1], 1.0) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<Fanourgakis::SP3Self> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2] * 2.0) );
                this->energyFunc.back()->set_geo(this->geo);

                Fanourgakis::R = args[0];

                printf("\tResetting box size in z to %lf\n", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
                this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2];
                break;

            case 13:
                printf("\nAdding Ewald with vacuum slabs\n");
                assert(args.size() == 6);

                printf("\tAdding vacuum slabs of thickness: %lf on each side of the box.\n", args[1]);
                this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 

                this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::Short> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::Long> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
                EwaldLike::alpha = args[5];
                printf("\tk-vectors: %d %d %d\n", (int) args[2], (int) args[3], (int) args[4]);
                break;

            case 14:
                printf("\nAdding Ewald slab correction\n");
                assert(args.size() == 0);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::SlabCorr> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                break;

            case 15:
                printf("\nAdding 2D Ewald potential\n");
                assert(args.size() == 5);
                EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
                EwaldLike::alpha = args[4];
                this->energyFunc.push_back( std::make_shared< Energy2D<EwaldLike::Ewald2D> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                printf("\tk-vectors: %d %d %d\n", (int) args[1], (int) args[2], (int) args[3]);
                break;

            case 16:
                printf("\nAdding Ewald slab correction2\n");
                assert(args.size() == 0);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::SlabCorr2> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                break;


            case 17:
                printf("\nAdding Ewald with vacuum slabs with plane-wise summation\n");
                assert(args.size() == 6);

                printf("\tAdding vacuum slabs of thickness: %lf on each side of the box.\n", args[1]);
                this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 

                this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::Short> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::Long2> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
                EwaldLike::alpha = args[5];
                printf("\tk-vectors: %d %d %d\n", (int) args[2], (int) args[3], (int) args[4]);
                break;

            case 18:
                printf("\nAdding Ewald slab correction3\n");
                assert(args.size() == 0);
                printf("\tBox size in z is: %lf\n", this->geo->d[2]);
                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::SlabCorr3> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                break;

            case 19:
                printf("\nAdding Repulsive Image Lennard-Jones\n");
                assert(args.size() == 1);
                this->energyFunc.push_back( std::make_shared< PairEnergyCOM<LJRep> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                break;

            case 20:
                printf("\nAdding Repulsive wall Lennard-Jones\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< ExternalEnergy<LJWallRep> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                break;

            case 21:
                printf("\nAdding Lennard-Jones\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< PairEnergyCOM<LJ> >(args[0], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                break;

            case 22:
                printf("\nAdding Repulsive wall Exponential\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< ExternalEnergy<ExpWallRep> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                break;
            case 23:
                printf("\nAdding Jan-potential to charges\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< ChargeWell<Jan> >(args[0], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(100.0);
                break;
            case 24:
                printf("\nAdding truncated and shifted LJ\n");
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< PairEnergyCOM<LJST> >(args[0], args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[1]);
                break;
            case 25:
                printf("\nAdding Ewald charge correction\n");
                assert(args.size() == 0);
                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::ChargeCorr> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                break;
            case 26:
                printf("\nAdding Ewald with vacuum slabs + correction terms (charged system)\n");
                assert(args.size() == 6);

                printf("\tAdding vacuum slabs of thickness: %lf on each side of the box.\n", args[1]);
                this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 

                this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::Short> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongChargedVacuum> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
                EwaldLike::alpha = args[5];
                printf("\tk-vectors: %d %d %d\n", (int) args[2], (int) args[3], (int) args[4]);
                break;
            case 27:
                printf("\nAdding Long range Ewald with vacuum slabs and excplicit wall charges\n");
                assert(args.size() == 6);

                printf("\tAdding vacuum slabs of thickness: %lf on each side of the box.\n", args[1]);
                this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
                this->_old->geo->d[2] = this->geo->d[2];
                this->_old->geo->dh[2] = this->geo->dh[2]; 

                this->energyFunc.push_back( std::make_shared< ExplicitWallChargeExtEnergy <EwaldLike::LongWithExplicitWallCharges, EwaldLike::Short> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                
                EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
                EwaldLike::alpha = args[5];
                printf("\tk-vectors: %d %d %d\n", (int) args[2], (int) args[3], (int) args[4]);
                break;
            case 28:
                printf("\nAdding Real part of Ewald with vacuum slabs\n");
                assert(args.size() == 1);

                this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::Short> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                break;
            default:
                printf("\nAdding Coulomb potential\n");
                this->energyFunc.push_back( std::make_shared< PairEnergy<Coulomb> >() );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);
                break;   
        }
        
    }

    void load_spline(std::vector<double> aKnots, std::vector<double> bKnots, std::vector<double >controlPoints){
        spline.load(aKnots, bKnots, controlPoints);
    }

    void load_cp(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<double> b_min, std::vector<double> b_max, std::vector<std::string> names){
        //assert correct sizes

        Eigen::Vector3d ae, be, qDisp;
        
        for(unsigned int i = 0; i < pos.size(); i++){
            ae << pos[i][0], pos[i][1], pos[i][2];
            be << com[i][0], com[i][1], com[i][2];
            qDisp = this->geo->displacement(ae, be);

            this->particles.add(com[i], pos[i], qDisp, r[i], rf[i], charges[i], b[i], b_min[i], b_max[i], names[i]);

        }
        if(!this->particles.setPModel || !this->particles.setNModel){
            printf("Cation or anion model not set!\n");
            exit(1);
        }
        printf("Loaded %u particles, %u cations and %u anions.\n", this->particles.tot, this->particles.cTot, this->particles.aTot);
    }

    void load_cp_old(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<std::string> names){
        //assert correct sizes

        Eigen::Vector3d ae, be, qDisp;
        double b_min, b_max;

        for(unsigned int i = 0; i < pos.size(); i++){
            ae << pos[i][0], pos[i][1], pos[i][2];
            be << com[i][0], com[i][1], com[i][2];
            qDisp = this->geo->displacement(ae, be);

            if(charges[i] < 0){
                b_min = 0.0;
                b_max = 0.0;
            }

            this->particles.add(com[i], pos[i], qDisp, r[i], rf[i], charges[i], b[i], b_min, b_max, names[i]);

        }
        if(!this->particles.setPModel || !this->particles.setNModel){
            printf("Cation or anion model not set!\n");
            exit(1);
        }
        printf("Loaded %u particles, %u cations and %u anions.\n", this->particles.tot, this->particles.cTot, this->particles.aTot);
    }
};