#define UNUSED(x) (void)(x)

#include "state.h"

namespace Simulator{


void State::advance(){
    this->step++;
}

void State::advanceMicro(){
    // Clear moved particles before next iteration
    this->movedParticles.clear();
    this->_old->movedParticles.clear();
}

void State::control(){
    #ifdef _DEBUG_
    Logger::Log<Logger::LogLevel::DEBUG>("Control (DEBUG)");
    #else
    Logger::Log("Control\n");
    #endif
    
    this->energy = 0.0;
    for(auto& e : this->energyFunc){
        if(this->step % 100 == 0){
            e->initialize(this->particles);
        }
        this->energy += e->all2all(this->particles);
    }

    // Check relative energy drift
    this->error = std::fabs((this->energy - this->cummulativeEnergy) / this->energy);
    if(this->error > 1e-10 || this->energy > 1e30){
        Logger::Log<Logger::LogLevel::FATAL>("Energy drift is too large: ", this->error, " (all2all: ", this->energy, 
                                                                    " cummulative: " ,this->cummulativeEnergy, ")");
        exit(0);
    } 

    // Debug checking
    #ifdef DEBUG
    if(this->particles.tot != _old->particles.tot){
        Logger::Log<Logger::LogLevel::FATAL>("_old state has ", _old->particles.tot, " particles and current state has", 
                                                                                            this->particles.tot);
        exit(0);
    }
    unsigned int cations = 0, anions = 0;
    for(const auto& p : this->particles){
        if(std::abs(this->geo->distance(p.com, p.pos) - p.b) > 1e-5){
            Logger::Log<Logger::LogLevel::FATAL>("|Pos - com| is not b! it is: ", 
                                                    this->geo->distance(p.com, p.pos), "and should be: ", p.b);
            exit(0);
        }
        if(std::abs(this->geo->distance(p.com, p.pos) - p.qDisp.norm()) > 1e-5){
            Logger::Log<Logger::LogLevel::FATAL>("|Pos - com| is not equal to |qDisp|!");
            exit(0);
        }
        if(p.b > p.b_max){
            Logger::Log<Logger::LogLevel::FATAL>("Ooops, b is larger than b_max!");
            exit(0);
        }
        if(p.pos != this->_old->particles[i]->pos){
            Logger::Log<Logger::LogLevel::FATAL>("current positions is not equal to old positions for particle ", i);
            std::cout << p.pos) << std::endl;
            std::cout << "\n" << this->_old->particles[i]->pos<< std::endl;
            exit(0);
        }
        if(p.com != this->_old->particles[i]->com){
            Logger::Log<Logger::LogLevel::FATAL>("current center of mass is not equal to old for particle ", i);
            std::cout << p.com << std::endl;
            std::cout << "\n" << this->_old->particles[p.index]->com << std::endl;
            exit(0);
        }
        if(p.index != i){
            Logger::Log<Logger::LogLevel::FATAL>("index is wrong in current for particle ", i, " it has index ", 
                                                                            p.index);
            exit(0);
        }
        if(this->_old->particles[i]->index != i){
            Logger::Log<Logger::LogLevel::FATAL>("index is wrong in in _old for particle ", i, " it has index ", 
                                                                        this->_old->particles[i]->index);
            exit(0);
        }
        (p.q > 0.0) ? cations++ : anions++;
    }

    if(cations != this->particles.cTot || anions != this->particles.aTot){
        Logger::Log<Logger::LogLevel::FATAL>("Wrong number of cations or anions!\n");
        exit(0);
    }

    if(this->particles.cTot + this->particles.aTot != this->particles.tot){
        Logger::Log<Logger::LogLevel::FATAL>("cTot + aTot != tot");
        exit(0);
    }

    if(!this->movedParticles.empty()){
        Logger::Log<Logger::LogLevel::FATAL>("Moved particles is not empty!");
        exit(0);
    }

    if(!this->_old->movedParticles.empty()){
        Logger::Log<Logger::LogLevel::FATAL>("_old Moved particles is not empty!");
        exit(0);
    }
    #endif
}

void State::finalize(std::string name){
    printf("\n");
    Logger::Log("Finalizing simulation: ", name.c_str());

    // Set up the _old system
    for(const auto& p : this->particles)
        _old->particles.add(p);

    //Calculate the initial energy of the system
    for(auto& e : this->energyFunc){
        e->initialize(particles);
        this->energy += e->all2all(this->particles);
    }

    this->cummulativeEnergy = this->energy;
    Logger::Log("\tEnergy of the first frame is: ", this->energy);
}

void State::reset_energy(){
    this->energy = 0.0;

    for(const auto& e : this->energyFunc)
        this->energy += e->all2all(this->particles);

    this->cummulativeEnergy = this->energy;
}

void State::save(){
    for(const auto i : this->movedParticles){
        // If a particle has been added, add it to the old state
        if(this->particles.tot > this->_old->particles.tot)
            this->_old->particles.add(this->particles[i]);

        else if(this->particles.tot == this->_old->particles.tot){
            *(this->_old->particles[i]) = *(this->particles[i]);
            //For SingleSwap move
            this->_old->particles.cTot = this->particles.cTot;
            this->_old->particles.aTot = this->particles.aTot;
        }
    }

    for(const auto i : this->_old->movedParticles){
        if(this->particles.tot < this->_old->particles.tot)
            this->_old->particles.remove(i);
    }

    this->cummulativeEnergy += this->dE;

    // If a volume move has been performed
    if(this->geo->volume != this->_old->geo->volume){
        //Update old geometry
        this->_old->geo->d      = this->geo->d;
        this->_old->geo->_d     = this->geo->_d;
        this->_old->geo->dh     = this->geo->dh;
        this->_old->geo->_dh    = this->geo->_dh;
        this->_old->geo->volume = this->geo->volume;
    }
}


void State::revert(){
    if(this->dE != std::numeric_limits<double>::infinity()){
        // Revert energy functors
        for(auto& e : this->energyFunc){
            if(this->geo->volume != this->_old->geo->volume){
                // Is update really needed here?
                e->update(this->_old->geo->d[0], this->_old->geo->d[1], this->_old->geo->d[2]);
                e->initialize(this->_old->particles);
            }
            else
                e->update( this->particles.get_subset(this->movedParticles), 
                                                    this->_old->particles.get_subset(this->_old->movedParticles) );
        }
    }

    //Set moved partiles in current state equal to previous state
    if(this->particles.tot == this->_old->particles.tot){
        for(const auto i : this->movedParticles){
            *(this->particles[i]) = *(this->_old->particles[i]);
            //For SingleSwap move
            this->particles.cTot = this->_old->particles.cTot;
            this->particles.aTot = this->_old->particles.aTot;
        }
    }
    // A particle has been added, remove it
    else if(this->particles.tot > this->_old->particles.tot){
        std::reverse(this->movedParticles.begin(), this->movedParticles.end());
        for(const auto i : this->movedParticles)
            this->particles.remove(i);
    }
    // A particle has been removed, add it back
    else if(this->particles.tot < this->_old->particles.tot){
        std::reverse(this->_old->movedParticles.begin(), this->_old->movedParticles.end());

        for(const auto i : this->_old->movedParticles)
            this->particles.add(this->_old->particles[i], i);
    }

    //Revert geometry
    this->geo->d      = this->_old->geo->d;
    this->geo->_d     = this->_old->geo->_d;
    this->geo->dh     = this->_old->geo->dh;
    this->geo->_dh    = this->_old->geo->_dh;
    this->geo->volume = this->_old->geo->volume;
}


//Get energy difference between this and old state
double State::get_energy_change(){
    auto E1 = 0.0, E2 = 0.0;

    for(const auto& p : this->movedParticles){
        // If the particle is outside of the box, or if there is hard sphere overlap
        if(!this->geo->is_inside(this->particles[p]) || this->overlap(p)){
            this->dE = std::numeric_limits<double>::infinity();
            return this->dE;
        }
    }

    for(auto& e : this->energyFunc){
        // Calculate the energy of the old state
        e->geo = this->_old->geo;
        E1 += (*e)( this->_old->movedParticles, this->_old->particles );

        // Calculate the energy of the new state
        e->geo = this->geo;
        if(this->geo->volume != this->_old->geo->volume){
            e->update(this->geo->d[0], this->geo->d[1], this->geo->d[2]);
            e->initialize(particles);
        }
        else
            e->update( this->_old->particles.get_subset(this->_old->movedParticles), 
                                                        this->particles.get_subset(this->movedParticles) );

        E2 += (*e)( this->movedParticles, this->particles );
    }

    this->dE = E2 - E1;
    return this->dE;
}


//Called after move - set movedParticles
void State::move_callback(std::vector< unsigned int > ps){  
    //If a particle is removed, this->movedparticles is empty. 
    //If particle is added this->_old->particles is empty
    if(this->particles.tot >= this->_old->particles.tot)
        std::for_each(std::begin(ps), std::end(ps), [this](int i){
                                                this->movedParticles.push_back(i); });
    
    std::copy_if(ps.begin(), ps.end(), std::back_inserter(this->_old->movedParticles), 
                                        [this](unsigned int i){ return i < this->_old->particles.tot; });

    for(auto p : this->movedParticles)
        geo->pbc(this->particles[p]);
}


void State::equilibrate(double step){
    printf("\n");
    Logger::Log("Equilibrating:");
    
    // Initial Check
    auto overlaps = this->get_overlaps();
    if(overlaps > 0){
        Logger::Log("\tRandomly placing particles");
        for(auto& p : this->particles){
            p->com = this->geo->random_pos(p->rf);
            p->pos = p->com + p->qDisp;
        }
    }

    Logger::Log("\tInitial overlaps: ", overlaps);

    unsigned long i = 0;
    //Move particles to prevent overlap
    if(overlaps > 0){
        Logger::Log("\tRemoving overlaps.");
        while(overlaps > 0){
            auto p = this->particles.random();
            Eigen::Vector3d oldCom = p->com;
            Eigen::Vector3d oldPos = p->pos;
            auto step_rand = Random::get_random() * step;
            p->translate(step_rand);
            //this->geo->pbc(p);
            if(!this->geo->is_inside(p) || this->overlap(p->index)){
                p->com = oldCom;
                p->pos = oldPos;
            }


            if(i % 50000 == 0){
                overlaps = this->get_overlaps();
                Logger::Log("\tOverlaps: ", overlaps, " iteration: ", i);
                fflush(stdout);
            }

            // This can take some time for dense systems so overflow is not impossible
            i++;
            if(i > 1E9)
                i = 0;
        }
    }
    else
        Logger::Log("\tNo overlaps to remove!");

    Logger::Log("\tEquilibration done");
}



bool State::overlap(std::size_t i) const{
    for(auto p : this->particles){
        if(p->index == i) continue;
        // Immediately return true if there is any overlap
        if(this->geo->distance(p->com, this->particles[i]->com) <= p->r + this->particles[i]->r)
            return true;
    }

    return false;
}

int State::get_overlaps() const{
    int count = 0;
    for(auto p : this->particles)
        (this->overlap(p->index)) ? count++ : 0;

    return count;
}

void State::set_geometry(int type, std::vector<double> args){
    this->_old = std::make_shared<State>();

    switch (type){
        default:
            Logger::Log("Creating Cuboid box");
            assert(args.size() == 3);
            this->geo = std::make_shared< Cuboid<true, true, true> >(args[0], args[1], args[2]);
            this->_old->geo = std::make_shared< Cuboid<true, true, true> >(args[0], args[1], args[2]);
            break;

        case 1:
            this->geo = std::make_shared<Sphere>();
            break;

        case 2:
            Logger::Log("Creating Cuboid-Image box");
            assert(args.size() == 3);
            this->geo = std::make_shared< CuboidImg<true, true, true> >(args[0], args[1], args[2]);
            this->_old->geo = std::make_shared< CuboidImg<true, true, true> >(args[0], args[1], args[2]);
            break;

        case 3:
            Logger::Log("Creating Cuboid-Image box with no PBC in z");
            assert(args.size() == 3);
            this->geo = std::make_shared< CuboidImg<true, true, false> >(args[0], args[1], args[2]);
            this->_old->geo = std::make_shared< CuboidImg<true, true, false> >(args[0], args[1], args[2]);
            break;

        case 4:
            Logger::Log("Creating Cuboid box with no PBC");
            assert(args.size() == 3);
            this->geo = std::make_shared< Cuboid<false, false, false> >(args[0], args[1], args[2]);
            this->_old->geo = std::make_shared< Cuboid<false, false, false> >(args[0], args[1], args[2]);
            break;

        case 5:
            Logger::Log("Creating Cuboid box with no PBC in z");
            assert(args.size() == 3);
            this->geo = std::make_shared< Cuboid<true, true, false> >(args[0], args[1], args[2]);
            this->_old->geo = std::make_shared< Cuboid<true, true, false> >(args[0], args[1], args[2]);
            break;
    }

}



void State::set_energy(int type, std::vector<double> args){
    printf("\n");
    switch (type){
        case 1:
            Logger::Log("\nAdding Ewald potential\n");
            if(bool(args[5])){
                printf("");
                assert(args.size() == 7);
            }
            else{
                assert(args.size() == 6);
            }
            
            this->energyFunc.push_back( std::make_shared< PairEnergy<Potentials::EwaldLike::Short> >() );
            //this->energyFunc.push_back( std::make_shared< PairEnergyWithRep<EwaldLike::Short> >(1) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::Long> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
            Potentials::EwaldLike::alpha = args[4];
            Potentials::EwaldLike::spherical = bool(args[5]);
            Potentials::EwaldLike::kMax = args[6];

            Logger::Log("\tReciprocal cutoff type: ", Potentials::EwaldLike::spherical ? "Spherical" : "Cubic");
            printf("\n");
            if(Potentials::EwaldLike::spherical)
                Logger::Log("\t\tMaximum number of k-vectors: ", Potentials::EwaldLike::kMax);
            Logger::Log("\tReciprocal cutoff: ", Potentials::EwaldLike::kMax);
            Logger::Log("\tk-vectors: ", (int) args[1], (int) args[2], (int) args[3]);
            break;

        case 2:
            Logger::Log("Adding Halfwald potential");
            assert(args.size() == 5);
            this->energyFunc.push_back( std::make_shared< ImgEnergy<Potentials::EwaldLike::Short> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::LongHW> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
            Potentials::EwaldLike::alpha = args[4];
            Logger::Log("\tk-vectors: ", (int) args[1], (int) args[2], (int) args[3]);
            break;
        
        case 3:
            Logger::Log("\nAdding HalfwaldIPBC potential");
            assert(args.size() == 5);
            this->energyFunc.push_back( std::make_shared< ImgEnergy<Potentials::EwaldLike::Short> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::LongHWIPBC> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
            Potentials::EwaldLike::alpha = args[4];
            break;

        case 4:
            Logger::Log("\nAdding Ellipsoidal Ewald");
            assert(args.size() == 5);
            this->energyFunc.push_back( std::make_shared< Ellipsoid<Potentials::BSpline2D> >(spline.aKnots, spline.bKnots, spline.controlPoints) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::LongEllipsoidal> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
            Potentials::EwaldLike::alpha = args[4];
            break;

        case 5:
            Logger::Log("\nAdding Minimum Image Halfwald");
            assert(args.size() == 3);
            this->energyFunc.push_back( std::make_shared< MIHalfwald<Potentials::Coulomb> >(args[1], args[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            // Set box size due to reflections, needed for distance PBC
            Logger::Log("\tResetting box size in z to: ", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
            this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 

            break;


        case 6:
            Logger::Log("Adding Truncated Ewald potential");
            assert(args.size() == 8);
            //this->energyFunc.push_back( std::make_shared< PairEnergy<EwaldLike::ShortTruncated> >() );
            this->energyFunc.push_back( std::make_shared< PairEnergyWithRep<Potentials::EwaldLike::ShortTruncated> >(1) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::LongTruncated> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
            Potentials::EwaldLike::alpha = args[4];

            Potentials::EwaldLike::R = args[5];
            Potentials::EwaldLike::kMax = args[6];
            Potentials::EwaldLike::spherical = bool(args[7]);
            Logger::Log("\tSpherical cutoff: ", Potentials::EwaldLike::spherical ? "true\n" : "false\n");
            Logger::Log("\tReciprocal cutoff: ", Potentials::EwaldLike::kMax);
            Potentials::EwaldLike::eta = Potentials::EwaldLike::R * 1.0 / (std::sqrt(2.0) * Potentials::EwaldLike::alpha);
            break;

        case 7:
            Logger::Log("Adding Halfwald with real replicates");
            assert(args.size() == 7);
            this->energyFunc.push_back( std::make_shared< MIHalfwald<Potentials::EwaldLike::Short> >(args[1], args[6]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            Logger::Log("\tResetting box size in z to: ", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
            this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::LongHW> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2] * 2.0) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
            Potentials::EwaldLike::alpha = args[5];

            break;

        case 8:
            Logger::Log("Adding harmonic well to charges");
            assert(args.size() == 1);
            this->energyFunc.push_back( std::make_shared< ChargeWell<Potentials::Harmonic> >(args[0], 0.0) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(100.0);
            break;

        case 9:
            Logger::Log("Adding Sture-potential to charges");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< ChargeWell<Potentials::Sture> >(args[0], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(100.0);
            break;

        case 10:
            Logger::Log("Adding FENE-potential to charges");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< ChargeWell<Potentials::FENE> >(args[0], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(100.0);
            break;

        case 11:
            Logger::Log("Adding FanourgakisSP2 with image charges");
            assert(args.size() == 2);
                                                                                        // kMax     eps
            this->energyFunc.push_back( std::make_shared< MIHalfwald<Potentials::Fanourgakis::SP2> >(args[1], 1.0) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::Fanourgakis::SP2Self> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2] * 2.0) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::Fanourgakis::R = args[0];

            Logger::Log("\tResetting box size in z to: ", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
            this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 
            break;

        case 12:
            Logger::Log("Adding FanourgakisSP3 with image charges");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< MIHalfwald<Potentials::Fanourgakis::SP3> >(args[1], 1.0) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::Fanourgakis::SP3Self> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2] * 2.0) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::Fanourgakis::R = args[0];

            Logger::Log("\tResetting box size in z to: ", (4.0 * args[1] + 2.0) * this->geo->_d[2]);
            this->geo->d[2] = (4.0 * args[1] + 2.0) * this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2];
            break;

        case 13:
            Logger::Log("Adding Ewald with vacuum slabs");
            assert(args.size() == 6);

            Logger::Log("\tAdding vacuum slabs of thickness: ", args[1], " on each side of the box.");
            this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 

            this->energyFunc.push_back( std::make_shared< PairEnergy<Potentials::EwaldLike::Short> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::Long> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
            Potentials::EwaldLike::alpha = args[5];
            Logger::Log("\tk-vectors: ", (int) args[2], (int) args[3], (int) args[4]);
            break;

        case 14:
            Logger::Log("Adding Ewald slab correction.");
            assert(args.size() == 0);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::SlabCorr> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            break;

        case 15:
            Logger::Log("Adding 2D Ewald potential");
            assert(args.size() == 5);
            Potentials::EwaldLike::set_km({ (int) args[1], (int) args[2], (int) args[3] });
            Potentials::EwaldLike::alpha = args[4];
            this->energyFunc.push_back( std::make_shared< Energy2D<Potentials::EwaldLike::Ewald2D> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            Logger::Log("\tk-vectors: ", (int) args[1], (int) args[2], (int) args[3]);
            break;

        case 16:
            Logger::Log("Adding Ewald slab correction2");
            assert(args.size() == 0);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::SlabCorr2> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            break;


        case 17:
            Logger::Log("Adding Ewald with vacuum slabs with plane-wise summation");
            assert(args.size() == 6);

            Logger::Log("\tAdding vacuum slabs of thickness: ", args[1], " on each side of the box.");
            this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 

            this->energyFunc.push_back( std::make_shared< PairEnergy<Potentials::EwaldLike::Short> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::Long2> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
            Potentials::EwaldLike::alpha = args[5];
            Logger::Log("\tk-vectors: ", (int) args[2], (int) args[3], (int) args[4]);
            break;

        case 18:
            Logger::Log("Adding Ewald slab correction3.");
            assert(args.size() == 0);
            Logger::Log("\tBox size in z is: ", this->geo->d[2]);
            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::SlabCorr3> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            break;

        case 19:
            Logger::Log("Adding Repulsive Lennard-Jones.");
            assert(args.size() == 1);
            this->energyFunc.push_back( std::make_shared< PairEnergyCOM<Potentials::LJRep> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            break;

        case 20:
            Logger::Log("Adding Repulsive wall Lennard-Jones");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< ExternalEnergy<Potentials::LJWallRep> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            break;

        case 21:
            Logger::Log("Adding Lennard-Jones");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< PairEnergyCOM<Potentials::LJ> >(args[0], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            break;

        case 22:
            Logger::Log("Adding Repulsive wall Exponential.");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< ExternalEnergy<Potentials::ExpWallRep> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            break;
        case 23:
            Logger::Log("Adding Jan-potential to charges.");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< ChargeWell<Potentials::Jan> >(args[0], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(100.0);
            break;
        case 24:
            Logger::Log("Adding truncated and shifted LJ.");
            assert(args.size() == 2);
            this->energyFunc.push_back( std::make_shared< PairEnergyCOM<Potentials::LJST> >(args[0], args[1]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[1]);
            break;
        case 25:
            Logger::Log("Adding Ewald charge correction.");
            assert(args.size() == 0);
            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::ChargeCorr> >(this->geo->_d[0], this->geo->_d[1], this->geo->_d[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            break;
        case 26:
            Logger::Log("Adding Ewald with vacuum slabs + correction terms (charged system).");
            assert(args.size() == 6);

            Logger::Log("\tAdding vacuum slabs of thickness: ", args[1], " on each side of the box.");
            this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 

            this->energyFunc.push_back( std::make_shared< PairEnergy<Potentials::EwaldLike::Short> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);

            this->energyFunc.push_back( std::make_shared< ExtEnergy<Potentials::EwaldLike::LongChargedVacuum> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);

            Potentials::EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
            Potentials::EwaldLike::alpha = args[5];
            Logger::Log("\tk-vectors: ", (int) args[2], (int) args[3], (int) args[4]);
            break;
        case 27:
            Logger::Log("Adding Long range Ewald with vacuum slabs and excplicit wall charges");
            assert(args.size() == 6);

            Logger::Log("\tAdding vacuum slabs of thickness: ", args[1], " on each side of the box.\n");
            this->geo->d[2] = 2.0 * args[1] + this->geo->_d[2];
            this->geo->dh[2] = 0.5 * this->geo->d[2]; 
            this->_old->geo->d[2] = this->geo->d[2];
            this->_old->geo->dh[2] = this->geo->dh[2]; 

            this->energyFunc.push_back( std::make_shared< ExplicitWallChargeExtEnergy <Potentials::EwaldLike::LongWithExplicitWallCharges, Potentials::EwaldLike::Short> >(this->geo->d[0], this->geo->d[1], this->geo->d[2]) );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            
            Potentials::EwaldLike::set_km({ (int) args[2], (int) args[3], (int) args[4] });
            Potentials::EwaldLike::alpha = args[5];
            Logger::Log("\tk-vectors: ", (int) args[2], (int) args[3], (int) args[4]);
            break;
        case 28:
            Logger::Log("Adding Real part of Ewald with vacuum slabs.");
            assert(args.size() == 1);

            this->energyFunc.push_back( std::make_shared< PairEnergy<Potentials::EwaldLike::Short> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            break;
        default:
            Logger::Log("Adding Coulomb potential.");
            this->energyFunc.push_back( std::make_shared< PairEnergy<Potentials::Coulomb> >() );
            this->energyFunc.back()->set_geo(this->geo);
            this->energyFunc.back()->set_cutoff(args[0]);
            break;   
    }
}

void State::load_spline(std::vector<double> aKnots, std::vector<double> bKnots, std::vector<double >controlPoints){
    spline.load(aKnots, bKnots, controlPoints);
}

void State::load_cp(std::vector< std::vector<double> > com, std::vector< std::vector<double> > pos, std::vector<double> charges, std::vector<double> r, std::vector<double> rf, std::vector<double> b, std::vector<double> b_min, std::vector<double> b_max, std::vector<std::string> names){
    //assert correct sizes
    printf("\n");
    Eigen::Vector3d ae, be, qDisp;
    
    for(unsigned int i = 0; i < pos.size(); i++){
        ae << pos[i][0], pos[i][1], pos[i][2];
        be << com[i][0], com[i][1], com[i][2];
        qDisp = this->geo->displacement(ae, be);

        this->particles.add(com[i], pos[i], qDisp, r[i], rf[i], charges[i], b[i], b_min[i], b_max[i], names[i]);

    }
    if(!this->particles.is_valid()){
        Logger::Log<Logger::LogLevel::FATAL>("Particles object is not valid!");
        exit(1);
    }
    Logger::Log("Loaded ", this->particles.tot, " particles, ", this->particles.cTot, " cations and ", this->particles.aTot, " anions.");
}


} // end of namespace Simulator