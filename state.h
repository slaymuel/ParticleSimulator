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

    std::shared_ptr<State> _old;
    SplineData spline;
    //IO io;
    
    public:

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
            this->energy += e->all2all(this->particles);;
        }

        this->error = std::fabs((this->energy - this->cummulativeEnergy) / this->energy);


        #ifdef DEBUG
        if(this->particles.tot != _old->particles.tot){
            printf("_old state has %i particles and current state has %i.\n", _old->particles.tot, this->particles.tot);
            exit(1);
        }

        for(unsigned int i = 0; i < this->particles.tot; i++){
            if(std::abs(this->geo->distance(this->particles.particles[i]->com, this->particles.particles[i]->pos) - 
                                                                this->particles.particles[i]->b) > 1e-5){
                printf("|Pos - com| is not b! it is: %.6lf and should be: %.6lf\n",
                                    this->geo->distance(this->particles.particles[i]->com, this->particles.particles[i]->pos), 
                                    this->particles.particles[i]->b);
                exit(1);
            }
            if(this->particles.particles[i]->pos != this->_old->particles.particles[i]->pos){
                printf("current positions is not equal to old positions for particle %i.\n", i);
                std::cout << this->particles.particles[i]->pos << std::endl;
                std::cout << "\n" << this->_old->particles.particles[i]->pos<< std::endl;
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
        }

        if(this->particles.cTot + this->particles.aTot != this->particles.tot){
            printf("Ctot + aTot != tot\n");
            exit(1);
        }
        #endif

        if(this->energy != 0 && this->cummulativeEnergy != 0){
            if(this->error > 1e-10 || this->energy > 1e30){
                printf("\n\nEnergy drift is too large: %.12lf (all2all: %lf, cummulative: %lf)\n\n", this->error, this->energy, this->cummulativeEnergy);
                exit(1);
            } 
        }
    }

    void finalize(std::string name){

        // Set up old system
        for(std::shared_ptr<Particle> p : this->particles.particles){
            _old->particles.add(p);
        }

        //Calculate the initial energy of the system
        for(auto e : this->energyFunc){
            e->initialize(particles);
            this->energy += e->all2all(this->particles);
        }
        this->cummulativeEnergy = this->energy;

        //io.open(name);
    }

    void reset_energy(){
        this->energy = 0.0;

        for(auto e : this->energyFunc){
            this->energy += e->all2all(this->particles);
        }

        this->cummulativeEnergy = this->energy;
    }

    void save(){
        //printf("Save:\n");
        //printf("current moved %lu\n", this->movedParticles.size());
        //printf("old moved %lu\n\n", this->_old->movedParticles.size());
        //printf("p1 %.8lf %.8lf %.8lf\n", this->particles[0]->com[0], this->particles[0]->com[1], this->particles[0]->com[2]);
        //printf("p1 old %.8lf %.8lf %.8lf\n", this->_old->particles[0]->com[0], this->_old->particles[0]->com[1], this->_old->particles[0]->com[2]);
        for(auto i : this->movedParticles){
            if(this->particles.tot > this->_old->particles.tot){
                //printf("\nSave: adding particle %i to old\n\n", i);
                this->_old->particles.add(this->particles.particles[i]);
                //this->_old->particles.add(this->particles.particles[i], i);
            }
            else if(this->particles.tot == this->_old->particles.tot){
                *(this->_old->particles.particles[i]) = *(this->particles.particles[i]);
                //For SingleSwap move
                this->_old->particles.cTot = this->particles.cTot;
                this->_old->particles.aTot = this->particles.aTot;
            }
        }
        //printf("p1 %.8lf %.8lf %.8lf\n", this->particles[0]->com[0], this->particles[0]->com[1], this->particles[0]->com[2]);
        //printf("p1 old %.8lf %.8lf %.8lf\n", this->_old->particles[0]->com[0], this->_old->particles[0]->com[1], this->_old->particles[0]->com[2]);
        for(auto i : this->_old->movedParticles){
            if(this->particles.tot < this->_old->particles.tot){
                //printf("\nSave: removing particle %i from old\n\n", i);
                this->_old->particles.remove(i);
            }
        }

        this->movedParticles.clear();
        this->_old->movedParticles.clear();
        this->cummulativeEnergy += this->dE;

        //printf("save b: old box: %lf %lf %lf\n", this->_old->geo->dh[0], this->_old->geo->dh[1], this->_old->geo->dh[2]);
        if(this->geo->volume != this->_old->geo->volume){
            //Update old geometry
            this->_old->geo->d = this->geo->d;
            this->_old->geo->_d = this->geo->_d;
            this->_old->geo->dh = this->geo->dh;
            this->_old->geo->_dh = this->geo->_dh;
            this->_old->geo->volume = this->geo->volume;
        }
        //printf("Save a: old box: %lf %lf %lf\n", this->_old->geo->dh[0], this->_old->geo->dh[1], this->_old->geo->dh[2]);
    }


    void revert(){
        //Set moved partiles in current state equal to previous state
        if(this->dE != std::numeric_limits<double>::infinity()){
            for(auto e : this->energyFunc){
                if(this->geo->volume != this->_old->geo->volume){
                    //printf("Volume move, reverting energy\n");
                    e->update(this->_old->geo->d[0], this->_old->geo->d[1], this->_old->geo->d[2]);
                    e->initialize(this->_old->particles);
                }
                else{
                    e->update( this->particles.get_subset(this->movedParticles), this->_old->particles.get_subset(this->_old->movedParticles) );
                }
            }
        }

        //printf("p1 %.8lf %.8lf %.8lf\n", this->particles[0]->com[0], this->particles[0]->com[1], this->particles[0]->com[2]);
        for(auto i : this->movedParticles){
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
        //printf("p1 %.8lf %.8lf %.8lf\n", this->particles[0]->com[0], this->particles[0]->com[1], this->particles[0]->com[2]);
        //   REARRANGE, MOVE CONDITION OUTSIDE OF LOOP, LESS GENERAL THOUGH.....
        for(auto i : this->_old->movedParticles){
            if(this->particles.tot < this->_old->particles.tot){
                //printf("\nRevert: adding back particle %i to current\n", i);
                this->particles.add(this->_old->particles.particles[i], i);
            }

        }

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

        for(auto p : this->movedParticles){
            if(!this->geo->is_inside(this->particles.particles[p]) || this->overlap(p)){
                this->dE = std::numeric_limits<double>::infinity();
                return this->dE;
            }
        }

        for(auto e : this->energyFunc){
            //stupid design
            e->geo = this->_old->geo;
            //auto start = std::chrono::steady_clock::now();

            E1 += (*e)( this->_old->movedParticles, this->_old->particles );
            //printf("b Energy: %lf\n", (*e)( this->_old->movedParticles, this->_old->particles ));
            //auto end = std::chrono::steady_clock::now();
            //std::cout << "Energy: " << (double) std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0 << "us\n\n";
            //start = std::chrono::steady_clock::now();
            e->geo = this->geo;
            if(this->geo->volume != this->_old->geo->volume){
                e->update(this->geo->d[0], this->geo->d[1], this->geo->d[2]);
                e->initialize(particles);
            }
            else{
                e->update( this->_old->particles.get_subset(this->_old->movedParticles), this->particles.get_subset(this->movedParticles) );
            }
            //printf("a Energy:%lf\n", (*e)( this->_old->movedParticles, this->_old->particles ));
            //end = std::chrono::steady_clock::now();
            //std::cout << "Update: " << (double) std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0 << "us " <<  "\n\n";
            E2 += (*e)( this->movedParticles, this->particles );
            //std::cout << this->_old->particles[this->_old->movedParticles[0]]->q << " " << this->particles[this->_old->movedParticles[0]]->q  << std::endl;
        }
        this->dE = E2 - E1;
        //printf("dE = %lf\n", this->dE);
        return this->dE;
    }


    //Called after move - set movedParticles
    void move_callback(std::vector< unsigned int > ps){   
        for(auto i : ps){
            geo->pbc(particles[i]);
        }

        //this->movedParticles.insert(std::end(movedParticles), std::begin(ps), std::end(ps));

        //If a particle is removed, this->movedparticles is empty. If particle is added this->_old->particles is empty
        if(this->particles.tot >= this->_old->particles.tot){
            std::for_each(std::begin(ps), std::end(ps), [this](int i){ 
                                                    this->movedParticles.push_back(i); });
        }

        std::copy_if(ps.begin(), ps.end(), std::back_inserter(this->_old->movedParticles), [this](int i){ return i < this->_old->particles.tot; });

        //printf("%lu particles moved, old %lu\n", this->movedParticles.size(),this->_old->movedParticles.size());
    }


    void equilibrate(double step){
        printf("\nEquilibrating:\n");
        Eigen::Vector3d v;
        
        // Initial Check
        int i = 0, overlaps = this->get_overlaps();
        //exit(0);
        if(overlaps > 0){
            for(unsigned int i = 0; i < this->particles.pTot; i++){
                this->particles.particles[i]->com = this->geo->random_pos(this->particles.particles[i]->rf);
                this->particles.particles[i]->pos = this->particles.particles[i]->com + this->particles.particles[i]->qDisp;
            }
        }

        printf("\tInitial overlaps: %i\n", overlaps);
        Eigen::Vector3d oldCom;
        Eigen::Vector3d oldPos;
        std::shared_ptr<Particle> p;

        //Move particles to prevent overlap
        while(overlaps > 0){
            p = this->particles.random();
            oldCom = p->com;
            oldPos = p->pos;
            p->translate(step);
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
        printf("\n\tEquilibration done\n\n");
    }



    bool overlap(std::size_t i){
        for(auto p : this->particles.particles){
            if(p->index == i) continue;

            if(this->geo->distance(p->com, this->particles.particles[i]->com) <= p->r + this->particles.particles[i]->r){
                //printf("%lf\n", this->geo->distance(p->com, this->particles.particles[i]->com));
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
                assert(args.size() == 2);
                this->energyFunc.push_back( std::make_shared< MIHalfwald<Coulomb> >(args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                // Set box size due to reflections
                printf("\tResetting box size in z to %lf\n", (4.0 * args[0] + 2.0) * this->geo->_d[2]);
                this->geo->d[2] = (4.0 * args[0] + 2.0) * this->geo->_d[2];
                this->geo->dh[2] = 0.5 * this->geo->d[2]; 
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
                assert(args.size() == 6);
                this->energyFunc.push_back( std::make_shared< MIHalfwald<EwaldLike::Short> >(args[1]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(args[0]);

                this->energyFunc.push_back( std::make_shared< ExtEnergy<EwaldLike::LongHW> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
                this->energyFunc.back()->set_geo(this->geo);

                EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
                EwaldLike::alpha = args[5];
                break;

            case 8:
                printf("\nAdding harmonic well to charges\n");
                assert(args.size() == 1);
                this->energyFunc.push_back( std::make_shared< ChargeWell<Harmonic> >(args[0]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(100.0);
                break;

            case 9:
                printf("\nAdding Sture-potential to charges\n");
                assert(args.size() == 1);
                this->energyFunc.push_back( std::make_shared< ChargeWell<Sture> >(args[0]) );
                this->energyFunc.back()->set_geo(this->geo);
                this->energyFunc.back()->set_cutoff(100.0);
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

/*
    void close(){
        io.close();
    }

    void to_xtc(int step, int time){
        io.to_xtc(particles, geo->d, step, time);
    }

    void to_gro(std::string fileName){
        io.to_gro(fileName, particles, geo->d);
    }
*/
};