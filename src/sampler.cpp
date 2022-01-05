#include "sampler.h"


namespace Simulator{


void Samplers::Density::sample(State& state){
    for(unsigned int i = 0; i < state.particles.tot; i++){
        //printf("%lu %i\n", this->density.size(), (int) (particles.particles[i]->pos[d] + this->dh));
        if(state.particles.particles[i]->q > 0){
            pDens.at( (unsigned int) ( (state.particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
        }

        else{
            nDens.at( (unsigned int) ( (state.particles[i]->pos[d] + this->dh) / this->binWidth ) )++;
        }
    }
    this->samples++;
}

void Samplers::Density::save(){
    std::ofstream f ("p" + dim + "_" + this->filename + ".txt");
    if (f.is_open())
    {
        for(unsigned int i = 0; i < this->pDens.size(); i++){
            f << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 -  this->dh<< " " <<  
                    (double) this->pDens[i] / (this->xb * this->yb * this->binWidth * this->samples) << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";
    
    std::ofstream fi ("n" + dim + "_" + this->filename + ".txt");
    if (fi.is_open())
    {
        for(unsigned int i = 0; i < this->nDens.size(); i++){
            fi << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 -  this->dh << " " <<  
                    (double) this->nDens[i] / (this->xb * this->yb * this->binWidth * this->samples) << "\n";
        }
        fi.close();
    }
    else std::cout << "Unable to open file";
}

void Samplers::Density::close(){};









void Samplers::Energy::sample(State &state){
    energies.push_back(state.cummulativeEnergy);
}

void Samplers::Energy::save(){
    std::ofstream f (this->filename);
    if (f.is_open()){
        for(const auto& e : energies){
            f << std::fixed << std::setprecision(15) << e << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";
}

void Samplers::Energy::close(){};











void Samplers::WidomHS::sample(State& state){
    Eigen::Vector3d com = state.geo->random_pos(2.5);
    Eigen::Vector3d qDisp;
    qDisp << 0.0, 0.0, 0.0;
    com[2] = (Random::get_random() * 0.2 - 0.1) * state.geo->_dh[2];
    //std::cout << com[0] << " " << com[1] << " " << com[2] << std::endl;
    state.particles.add(com, com, qDisp, 2.5, state.particles.pModel.rf, state.particles.pModel.q, state.particles.pModel.b, 0.0, 0.0, "WIDOM_PARTICLE");

    if(!state.overlap(state.particles.tot - 1)){
        this->cp += 1.0;
    }
    state.particles.remove(state.particles.tot - 1);
    this->samples++;
    //printf("HS-CP: %lf\n", -std::log(this->cp / this->samples));
    //printf("samples: %d\n", this->samples);
    //printf("cp: %lf\n\n", this->cp);
}

void Samplers::WidomHS::save(){
    std::ofstream f (this->filename);
    if (f.is_open()){
        f << std::fixed << std::setprecision(10) << "Hard-sphere chemical potential: " << -std::log(this->cp / this->samples)  << "\n";
        f.close();
    }
    else std::cout << "Unable to open file";
}

void Samplers::WidomHS::close(){};






void Samplers::QDist::sample(State& state){
    for(unsigned int i = 0; i < state.particles.tot; i++){
        int index = (int)((state.geo->distance(state.particles.particles[i]->pos, state.particles.particles[i]->com)) / this->binWidth);
        if(state.particles.particles[i]->q > 0.0){
            
            pqDist[index]++;
        }
        else{
            nqDist[index]++; 
        }
    }
    this->samples++;
}

void Samplers::QDist::save(){
    std::ofstream pf ("pqDist_" + this->filename + ".txt");
    if (pf.is_open())
    {
        for(unsigned int i = 0; i < this->pqDist.size(); i++){
            pf << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                    (double) this->pqDist[i] / this->samples << "\n";
        }
        pf.close();
    }
    else std::cout << "Unable to open file";

    std::ofstream nf ("nqDist_" + this->filename + ".txt");
    if (nf.is_open())
    {
        for(unsigned int i = 0; i < this->nqDist.size(); i++){
            nf << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                    (double) this->nqDist[i] / this->samples << "\n";
        }
        nf.close();
    }
    else std::cout << "Unable to open file";
}

void Samplers::QDist::close(){};









void Samplers::XDR::save(){};

void Samplers::XDR::close(){
    xdrfile_close(xdf);
}

void Samplers::XDR::sample(State& state){
    matrix box;
    box[0][0] = state.geo->_d[0];
    box[0][1] = 0.0;
    box[0][2] = 0.0;
    box[1][0] = state.geo->_d[1];
    box[1][1] = 0.0;
    box[1][2] = 0.0;
    box[2][0] = state.geo->_d[2];
    box[2][1] = 0.0;
    box[2][2] = 0.0;

    
    if (xdf != nullptr) {
        rvec *ps = new rvec[state.particles.tot];
        //size_t N = 0;

        for (unsigned int i = 0; i < state.particles.tot; i++) {
            ps[i][0] = state.particles[i]->pos[0] * 0.1 + state.geo->_d[0] * 0.5;
            ps[i][1] = state.particles[i]->pos[1] * 0.1 + state.geo->_d[1] * 0.5;
            ps[i][2] = state.particles[i]->pos[2] * 0.1 + state.geo->_d[2] * 0.5; 
        }

        write_xtc(xdf, state.particles.tot, state.step, state.step, box, ps, 1000);

        delete[] ps;
    }
    else{
        printf("Could not open xtc file!\n");
        exit(0);
    }
}







void Samplers::NumIons::sample(State& state){
    this->pNum.push_back(state.particles.cTot);
    this->nNum.push_back(state.particles.aTot);
}

void Samplers::NumIons::save(){
    std::ofstream f ("pNum_" + this->filename + ".txt");
    if (f.is_open())
    {
        for(auto p : this->pNum){
            f << p << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";
    
    std::ofstream fi ("nNum_" + this->filename + ".txt");
    if (fi.is_open())
    {
        for(auto n : this->nNum){
            fi << n << "\n";
        }
        fi.close();
    }
    else std::cout << "Unable to open file";
}

void Samplers::NumIons::close(){};






void Samplers::Pressure::sample(State& state){
    Eigen::Vector3d force;
    Eigen::Vector3d disp;
    Eigen::Vector3d temp;
    double press = 0.0;
    /*for(int i = 0; i < state.particles.tot; i++){
        force = Eigen::Vector3d::Zero();
        for(auto& e : state.energyFunc){
            force += (*e).force(state.particles[i], state.particles);
        }
        std::cout << force << std::endl;
    }*/

    #pragma omp parallel for private(force, disp, temp) reduction(+:this->pressureT, press)
    for(int i = 0; i < state.particles.tot; i++){
        for(int j = i + 1; j < state.particles.tot; j++){
            force = Eigen::Vector3d::Zero();
            for(const auto& e : state.energyFunc){
                force += (*e).force(state.particles[i], state.particles[j]);
            }
            //printf("Force: ");
            //std::cout << force << std::endl;
            //std::cout << force.norm() / constants::lB << std::endl;
            disp = state.geo->displacement(state.particles[i]->pos, state.particles[j]->pos);
            this->pressureT(0, 0) += disp[0] * force[0];   this->pressureT(0, 1) += disp[0] * force[1];   this->pressureT(0, 2) += disp[0] * force[2];
            this->pressureT(1, 0) += disp[1] * force[0];   this->pressureT(1, 1) += disp[1] * force[1];   this->pressureT(1, 2) += disp[1] * force[2];
            this->pressureT(2, 0) += disp[2] * force[0];   this->pressureT(2, 1) += disp[2] * force[1];   this->pressureT(2, 2) += disp[2] * force[2];
            press += disp[2] * force[2];

            

            
            temp = state.particles[j]->pos;
            temp[2] = math::sgn(state.particles[j]->pos[2]) * lZ - state.particles[j]->pos[2]; 
            
            force = Eigen::Vector3d::Zero();
            for(const auto& e : state.energyFunc){
                force += (*e).force(state.particles[i]->pos, temp, state.particles[i]->q, -state.particles[j]->q);
            }
            disp = state.geo->displacement(state.particles[i]->pos, temp);
            this->pressureT(0, 0) += disp[0] * force[0];   this->pressureT(0, 1) += disp[0] * force[1];   this->pressureT(0, 2) += disp[0] * force[2];
            this->pressureT(1, 0) += disp[1] * force[0];   this->pressureT(1, 1) += disp[1] * force[1];   this->pressureT(1, 2) += disp[1] * force[2];
            this->pressureT(2, 0) += disp[2] * force[0];   this->pressureT(2, 1) += disp[2] * force[1];   this->pressureT(2, 2) += disp[2] * force[2];
            press += disp[2] * force[2];
            
        }
    }
    pressures[counter] = press / state.geo->volume;
    this->idP += state.particles.tot / state.geo->volume;
    this->samples++;
    counter++;
}

void Samplers::Pressure::save(){
    printf("Ideal pressure: %lf\n", idP / this->samples);
    std::cout << "Pressure tensor:" << std::endl;
    std::cout <<  this->pressureT * 1.0 / this->V * 1.0 / this->samples << std::endl;
    std::cout << "Pressure in z: " << idP / this->samples + this->pressureT(2, 2) * 1.0 / this->V * 1.0 / this->samples << std::endl;

    std::ofstream f (this->filename + ".txt", std::ios::app);
    if (f.is_open())
    {
        for(int i = 0; i < counter; i++){
            f << std::fixed << std::setprecision(10) << this->pressures[i] << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";
    counter = 0;
}

void Samplers::Pressure::close(){};








void Samplers::PressureV::sample(State& state){
    bool fail = false;
    double e1 = 0.0, e2 = 0.0;
    oldd = state.geo->d[2];
    old_d = state.geo->_d[2];

    for(const auto& e : state.energyFunc){
        e1 += (*e).all2all(state.particles);
    }

    state.geo->_d[2] += this->dL;
    state.geo->_dh[2] = state.geo->_d[2] / 2.0;
    state.geo->d[2] += this->dL;
    state.geo->dh[2] = state.geo->d[2] / 2.0;
    state.geo->volume = this->newV;

    for(unsigned int i = 0; i < state.particles.tot; i++){
        state.particles[i]->com[2] *= state.geo->d[2] / oldd;
        state.particles[i]->pos = state.particles[i]->com + state.particles[i]->qDisp;
        state.geo->pbc(state.particles[i]);
        if(state.overlap(i)){
            fail = true;
        }
    }

    for(auto& e : state.energyFunc){
        //initialize reciprocal vectors again
        (*e).update(state.geo->d[0], state.geo->d[1], state.geo->d[2]);
        (*e).initialize(state.particles);
        e2 += (*e).all2all(state.particles);
    }

    for(unsigned int i = 0; i < state.particles.tot; i++){
        state.particles[i]->com[2] = state._old->particles[i]->com[2];
        state.particles[i]->pos[2] = state._old->particles[i]->pos[2];
    }

    state.geo->_d[2] = old_d;
    state.geo->_dh[2] = old_d / 2.0;
    state.geo->d[2] = oldd;
    state.geo->dh[2] = oldd / 2.0;
    state.geo->volume = oldV;

    for(auto& e : state.energyFunc){
        (*e).update(state.geo->d[0], state.geo->d[1], state.geo->d[2]);
        (*e).initialize(state.particles);
    }

    if(!fail){
        this->av += std::pow(1.0 + this->dV / state.geo->volume, state.particles.tot) * std::exp(-(e2 - e1));
        this->pressures[counter] = 1.0 / this->dV * std::log(std::pow(1.0 + this->dV / state.geo->volume, state.particles.tot) * std::exp(-(e2 - e1)));
        this->samples++;
        counter++;
    }
}

void Samplers::PressureV::save(){
    printf("PressureV: %lf\n", 1.0 / this->dV * std::log(this->av / this->samples));
    std::ofstream f (this->filename + ".txt", std::ios::app);
    if (f.is_open())
    {
        for(int i = 0; i < counter; i++){
            f << std::fixed << std::setprecision(10) << this->pressures[i] << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";

    counter = 0;
}

void Samplers::PressureV::close(){};









void Samplers::ForcePressure::sample(State& state){
    Eigen::Vector3d force;
    Eigen::Vector3d disp;
    Eigen::Vector3d temp;

    rightToLeft = 0.0;
    leftToRight = 0.0;
    
    //#pragma omp parallel for private(force, disp) reduction(+:leftToRight, +:rightToLeft)
    for(int i = 0; i < state.particles.tot; i++){
        for(int j = 0; j < state.particles.tot; j++){
            if(i==j)
                continue;

            int counter = 0;
            for(const auto& e : state.energyFunc){
                if(state.particles[i]->pos[2] < 0.0 && state.particles[j]->pos[2] >= 0.0){
                    ////Force on particle i from j
                    force = (*e).force(state.particles[i], state.particles[j]);

                    //temp is position of Image charge
                    temp = state.particles[j]->pos;
                    //lZ is length of C cell in z
                    temp[2] = math::sgn(temp[2]) * lZ - temp[2]; 

                    if(counter < 2){
                        force += (*e).force(state.particles[i]->pos, temp, state.particles[i]->q, -state.particles[j]->q);
                        rightToLeft += force[2];
                    }
                }

                else if(state.particles[j]->pos[2] < 0.0 && state.particles[i]->pos[2] >= 0.0){
                    //Force on particle j from i
                    force = (*e).force(state.particles[j], state.particles[i]);
                    temp = state.particles[i]->pos;
                    temp[2] = math::sgn(temp[2]) * lZ - temp[2]; 

                    if(counter < 2){
                        force += (*e).force(state.particles[j]->pos, temp, state.particles[j]->q, -state.particles[i]->q);
                        leftToRight += force[2];
                    }
                }
                counter++;
            }
        }
    }
    this->leftForce[this->counter] = rightToLeft;
    this->rightForce[this->counter] = leftToRight;
    this->idP += state.particles.tot / state.geo->volume;
    this->counter++;
    this->samples++;
}

void Samplers::ForcePressure::save(){
    printf("Ideal pressure: %lf\n", idP / this->samples);
    std::cout << "RightToLeft: " << this->leftForce[this->counter - 1] / (80.0 * 80.0) << std::endl;
    std::cout << "LeftToRight: " << this->rightForce[this->counter - 1] / (80.0 * 80.0) << std::endl;
    /*std::cout << "Right to left pressure: " << rightToLeft / rightToLeftSamples << std::endl;
    std::cout << "Right to left: " << rightToLeft / rightToLeftSamples << std::endl;
    std::cout << "Left to right: " << leftToRight / leftToRightSamples << std::endl;*/

    std::ofstream f (this->filename + ".txt", std::ios::app);
    if (f.is_open())
    {
        for(unsigned int i = 0; i < this->leftForce.size(); i++){
            f << std::fixed << std::setprecision(10) << this->samples - this->leftForce.size() + i << "\t" << leftForce[i] << "\t" << rightForce[i]  << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";

    this->counter = 0;
}

void Samplers::ForcePressure::close(){};










void Samplers::Force::sample(State& state){

    //#pragma omp parallel for private(force, disp) reduction(+:leftToRight, +:rightToLeft)
    for(int i = 0; i < state.particles.tot; i++){
        force = Eigen::Vector3d::Zero();
        for(int j = 0; j < state.particles.tot; j++){
            if(i==j)
                continue;

            for(const auto& e : state.energyFunc){
                    //Force that points towards particle i
                force += (*e).force(state.particles[i], state.particles[j]);
            }
        }
        std::cout << "Force is: " << force / constants::lB << std::endl;
    }
}

void Samplers::Force::save(){}
void Samplers::Force::close(){};













void Samplers::CliffPressure::sample(State& state){
    double rPTemp = 0.0;
    double lPTemp = 0.0;
    // Move particles
    for(int i = 0; i < state.particles.tot; i++){
        if(state.particles[i]->com[2] > 0.0){
            state.particles[i]->com[2] -= this->dL;
            state.particles[i]->pos = state.particles[i]->com + state.particles[i]->qDisp;
        }       
    }
    for(auto& e : state.energyFunc){
        //(*e).update(state.geo->d[0], state.geo->d[1], state.geo->d[2]);
        (*e).initialize(state.particles);
        rPTemp += (*e).all2all(state.particles);
    }
    for(int i = 0; i < state.particles.tot; i++){
        if(state.particles[i]->com[2] > 0.0){
            state.particles[i]->com[2] = state._old->particles[i]->com[2];
            state.particles[i]->pos = state._old->particles[i]->pos;
        }
    }
    for(auto& e : state.energyFunc){
        //(*e).update(state.geo->d[0], state.geo->d[1], state.geo->_d[2]);
        (*e).initialize(state.particles);
    }

    //Move right wall
    double oldd = state.geo->d[2];
    double old_d = state.geo->_d[2];
    state.geo->_d[2] += this->dL;
    state.geo->_dh[2] = 0.5*state.geo->_d[2];
    state.geo->d[2] = 2.0 * state.geo->_d[2];
    state.geo->dh[2] = 0.5*state.geo->d[2];
    for(int i = 0; i < state.particles.tot; i++){
        state.particles[i]->com[2] -= 0.5 * dL;
        state.particles[i]->pos = state.particles[i]->com + state.particles[i]->qDisp;
    }

    for(auto& e : state.energyFunc){
        (*e).update(state.geo->d[0], state.geo->d[1], state.geo->d[2]);
        (*e).initialize(state.particles);
        lPTemp += (*e).all2all(state.particles);
    }
    for(int i = 0; i < state.particles.tot; i++){
        state.particles[i]->com[2] = state._old->particles[i]->com[2];
        state.particles[i]->pos[2] = state._old->particles[i]->pos[2];
    }
    state.geo->_d[2] = old_d;
    state.geo->_dh[2] = 0.5*state.geo->_d[2];
    state.geo->d[2] = oldd;
    state.geo->dh[2] = 0.5*state.geo->d[2];
    for(auto& e : state.energyFunc){
        (*e).update(state.geo->d[0], state.geo->d[1], state.geo->d[2]);
        (*e).initialize(state.particles);
    }

    this->samples++;
    this->rP += rPTemp;
    this->lP += lPTemp;  
    this->pressures[counter] = (rPTemp - lPTemp) / (this->dL * this->area);
    counter++;
}

void Samplers::CliffPressure::save(){
    std::cout << "rP: " << this->rP / (this->samples * this->dL * this->area) << std::endl;
    std::cout << "lP: " << this->lP / (this->samples * this->dL * this->area) << std::endl;
    std::cout << "Cliff pressure: " << (this->rP - this->lP) / (this->samples * this->dL * this->area) << std::endl;

    std::ofstream f (this->filename + ".txt", std::ios::app);
    if (f.is_open())
    {
        for(unsigned int i = 0; i < this->pressures.size(); i++){
            f << std::fixed << std::setprecision(10) << this->samples - this->pressures.size() + i << "\t" << pressures[i] << "\n";
        }
        f.close();
    }
    else std::cout << "Unable to open file";

    this->counter = 0;
}

void Samplers::CliffPressure::close(){};









void Samplers::ModifiedWidom::sample(State& state){
    double elDE = 0.0, startExt;
    Eigen::Vector3d com = state.geo->random_pos(state.particles.pModel.rf);
    Eigen::Vector3d qDisp;
    qDisp << 0.0, 0.0, 0.0;

    //startExt = -state.energyFunc[0]->i2all(state.particles[state.particles.tot - 1], state.particles);
    //startExt -= state.energyFunc[2]->i2all(state.particles[state.particles.tot - 1], state.particles);

    std::vector< unsigned int > pIndice = {};
    pIndice.push_back(state.particles.tot);
    startExt = -(*state.energyFunc[0])(pIndice, state.particles);
    startExt -= (*state.energyFunc[2])(pIndice, state.particles);

    //Add cation
    state.particles.add(com, com, qDisp, state.particles.pModel.r, state.particles.pModel.rf, state.particles.pModel.q, 0.0, 0.0, 0.0, "WIDOM_PARTICLE");


    /////////////////////////////////////////////////////////////// Energy of inserting a full cation ///////////////////////////////////////////////////////////////////////////
    double sElDE = startExt;
    double sSrDE = 0.0;
    int k = 0;
    for(auto& e : state.energyFunc){
        std::vector< std::shared_ptr<Particle> > empty = {};
        e->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
        if(k < elecStart){
            //sElDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
            sElDE += (*e)(pIndice, state.particles);
            //printf("%i %lf\n", k, elDE);
        }
        else{
            sSrDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
            //printf("srDe %lf\n", srDE);
        }
        e->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
        k++;
    }    
    //sElDE *= constants::lB;
    ////////////////////////////////////////////////////////////// Integrate charge //////////////////////////////////////////////////////////////////////////////////////////////
    for(int scale = 1; scale <= 10; scale++){
        elDE = startExt;
        state.particles.particles[state.particles.tot - 1]->q = ((double) (scale - 1) * 0.1 + 0.1) * state.particles.pModel.q;
        
        std::vector< std::shared_ptr<Particle> > empty = {};
        int k = 0;
        for(auto& e : state.energyFunc){
            if(k < elecStart){
                e->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
                //elDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
                elDE += (*e)(pIndice, state.particles);
                e->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
            }
            k++;
        }    
        //elDE *= constants::lB;
        //nomP[scale] += 1.0 / ((double) (scale - 1) * 0.1 + 0.1) * elDE * std::exp(-(sSrDE + elDE));
        nomP[scale] += sElDE * std::exp(-(sSrDE + elDE));
        //std::cout << 1.0 / ((double) (scale - 1) * 0.1 + 0.1) * elDE << " " << sElDE << std::endl;
        denomP[scale] += std::exp(-(sSrDE + elDE));
    }

    nomP[0] += sElDE * std::exp(-sSrDE);
    denomP[0] += std::exp(-sSrDE);
    //if(!state.overlap(state.particles.tot - 1)){
        facP += std::exp(-sSrDE);
    //}
    pSamples++;

    state.particles.remove(state.particles.tot - 1);

    /////////////////////////////////////////////////////////////////////////////// Anion /////////////////////////////////////////////////////////////////////////////////////
    com = state.geo->random_pos(state.particles.nModel.rf);
    sElDE = startExt;
    sSrDE = 0.0;
    k = 0;
    state.particles.add(com, com, qDisp, state.particles.nModel.r, state.particles.nModel.rf, state.particles.nModel.q, 0.0, 0.0, 0.0, "WIDOM_PARTICLE");
    /////////////////////////////////////////////////////////////// Calculate full addition energy /////////////////////////////////////////////////////////////////////////////
    for(auto& e : state.energyFunc){
        std::vector< std::shared_ptr<Particle> > empty = {};
        e->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
        if(k < elecStart){
            sElDE += (*e)(pIndice, state.particles);
            //sElDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
        }
        else{
            sSrDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
            //printf("srDe %lf\n", srDE);
        }
        e->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
        k++;
    }   
    //sElDE *= constants::lB;
    ////////////////////////////////////////////////////////////// Integrate charge //////////////////////////////////////////////////////////////////////////////////////////////
    for(int scale = 1; scale <= 10; scale++){
        elDE = startExt;
        state.particles.particles[state.particles.tot - 1]->q = ((double) (scale - 1) * 0.1 + 0.1) * state.particles.nModel.q;
        std::vector< std::shared_ptr<Particle> > empty = {};
        int k = 0;
        for(auto& e : state.energyFunc){
            if(k < elecStart){
                e->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
                elDE += (*e)(pIndice, state.particles);
                //elDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
                e->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
            }
            k++;
        }  
        //elDE *= constants::lB;
        nomN[scale] += sElDE * std::exp(-(sSrDE + elDE));
        denomN[scale] += std::exp(-(sSrDE + elDE));  
    }
    nomN[0] += sElDE * std::exp(-sSrDE);
    denomN[0] += std::exp(-sSrDE);

    //if(!state.overlap(state.particles.tot - 1)){
        facN += std::exp(-sSrDE);
    //}
    nSamples++;
    
    state.particles.remove(state.particles.tot - 1);

    this->samples++;
}

void Samplers::ModifiedWidom::save(){
    double integralP = 0.0;
    for(int i = 0; i <= 10; i++){
        integralP += (nomP[i] / this->pSamples) / (denomP[i] / this->pSamples) * 1.0 / 11.0;
    }
    double integralN = 0.0;
    for(int i = 0; i <= 10; i++){
        integralN += (nomN[i] / this->nSamples) / (denomN[i] / this->nSamples) * 1.0 / 11.0;
    }


    double sumP = 0.5 * (nomP[0] / this->pSamples) / (denomP[0] / this->pSamples) + 0.5 * (nomP[10] / this->pSamples) / (denomP[10] / this->pSamples);
    double sumN = 0.5 * (nomN[0] / this->nSamples) / (denomN[0] / this->nSamples) + 0.5 * (nomN[10] / this->nSamples) / (denomN[10] / this->nSamples);

    printf("MW + sr: %lf el %lf tot %lf\n", -std::log(this->facP / this->pSamples), integralP, -std::log(this->facP / this->pSamples) + integralP);
    printf("MW - sr: %lf el %lf tot %lf\n", -std::log(this->facN / this->nSamples), integralN, -std::log(this->facN / this->nSamples) + integralN);
    printf("MWSum + el: %lf tot %lf\n", sumP, -std::log(this->facP / this->pSamples) + sumP);
    printf("MWSum - el: %lf tot %lf\n", sumN, -std::log(this->facN / this->nSamples) + sumN);
}

void Samplers::ModifiedWidom::close(){};













void Samplers::ModifiedWidomCoulomb::sample(State& state){
    e->set_geo(state.geo);
    double elDE = 0.0;
    Eigen::Vector3d com = state.geo->random_pos(state.particles.pModel.rf);
    Eigen::Vector3d qDisp;
    qDisp << 0.0, 0.0, 0.0;

    //Add cation
    state.particles.add(com, com, qDisp, state.particles.pModel.r, state.particles.pModel.rf, state.particles.pModel.q, 0.0, 0.0, 0.0, "WIDOM_PARTICLE");
    /////////////////////////////////////////////////////////////// Energy of inserting a full cation ///////////////////////////////////////////////////////////////////////////
    double sElDE = 0.0;
    double sSrDE = 0.0;
    int k = 0;

    for(int l = 0; l < state.particles.tot - 1; l++){
        state.particles.particles[l]->q *= (1.0 - state.particles.particles[l]->q / (state.particles.tot - 1));
    }
    sElDE += e->i2all(state.particles[state.particles.tot - 1], state.particles);
    for(int l = 0; l < state.particles.tot - 1; l++){
        state.particles.particles[l]->q = state.particles.particles[l]->q < 0.0 ? -1.0 : 1.0;
    }
    for(auto& e1 : state.energyFunc){

        std::vector< std::shared_ptr<Particle> > empty = {};
        e1->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
        if(k >= elecStart){
            sSrDE += e1->i2all(state.particles[state.particles.tot - 1], state.particles);
            //printf("srDe %lf\n", srDE);
        }
        empty = {};
        e1->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
        k++;
    }    
    sElDE *= constants::lB;
    ////////////////////////////////////////////////////////////// Integrate charge //////////////////////////////////////////////////////////////////////////////////////////////
    for(int scale = 1; scale <= 10; scale++){
        double lambda = ((double) (scale - 1.0) * 0.1 + 0.1);
        state.particles.particles[state.particles.tot - 1]->q = lambda * state.particles.pModel.q;
        for(int l = 0; l < state.particles.tot - 1; l++){
            state.particles.particles[l]->q *= (1.0 - state.particles.particles[l]->q * lambda / (state.particles.tot - 1));
        }

        std::vector< std::shared_ptr<Particle> > empty = {};
        e->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
        elDE = e->i2all(state.particles[state.particles.tot - 1], state.particles);
        empty = {};
        e->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
        elDE *= constants::lB;

        //nomP[scale] += 1.0 / ((double) (scale - 1) * 0.1 + 0.1) * elDE * std::exp(-(sSrDE + elDE));
        nomP[scale] += sElDE * std::exp(-(sSrDE + elDE));
        //std::cout << 1.0 / ((double) (scale - 1) * 0.1 + 0.1) * elDE << " " << sElDE << std::endl;
        denomP[scale] += std::exp(-(sSrDE + elDE));

        for(int l = 0; l < state.particles.tot - 1; l++){
            state.particles.particles[l]->q = state.particles.particles[l]->q < 0.0 ? -1.0 : 1.0;
        }
    }
    nomP[0] += elDE * std::exp(-sSrDE);
    denomP[0] += std::exp(-sSrDE);
    //if(!state.overlap(state.particles.tot - 1)){
        facP += std::exp(-sSrDE);
    //}
    pSamples++;

    state.particles.remove(state.particles.tot - 1);

    /////////////////////////////////////////////////////////////////////////////// Anion /////////////////////////////////////////////////////////////////////////////////////
    com = state.geo->random_pos(state.particles.nModel.rf);
    sElDE = 0.0;
    sSrDE = 0.0;
    k = 0;
    state.particles.add(com, com, qDisp, state.particles.nModel.r, state.particles.nModel.rf, state.particles.nModel.q, 0.0, 0.0, 0.0, "WIDOM_PARTICLE");
    /////////////////////////////////////////////////////////////// Calculate full addition energy /////////////////////////////////////////////////////////////////////////////
    for(int l = 0; l < state.particles.tot - 1; l++){
        state.particles.particles[l]->q *= (1.0 + state.particles.particles[l]->q / (state.particles.tot - 1));
    }
    sElDE = e->i2all(state.particles[state.particles.tot - 1], state.particles);
    for(int l = 0; l < state.particles.tot - 1; l++){
        state.particles.particles[l]->q = state.particles.particles[l]->q < 0.0 ? -1.0 : 1.0;
    }
    for(auto& e1 : state.energyFunc){
        std::vector< std::shared_ptr<Particle> > empty = {};
        e1->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
        if(k >= elecStart){
            sSrDE += e1->i2all(state.particles[state.particles.tot - 1], state.particles);
            //printf("srDe %lf\n", srDE);
        }
        empty = {};
        e1->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );
        k++;
    }   
    sElDE *= constants::lB;
    ////////////////////////////////////////////////////////////// Integrate charge //////////////////////////////////////////////////////////////////////////////////////////////
    for(int scale = 1; scale <= 10; scale++){
        double lambda = ((double) (scale - 1) * 0.1 + 0.1);
        state.particles.particles[state.particles.tot - 1]->q = lambda * state.particles.nModel.q;
        for(int l = 0; l < state.particles.tot - 1; l++){
            state.particles.particles[l]->q *= (1.0 + state.particles.particles[l]->q * lambda / (state.particles.tot - 1));
        }

        std::vector< std::shared_ptr<Particle> > empty = {};
        e->update( std::move(empty), state.particles.get_subset(state.particles.tot - 1) );
        elDE = e->i2all(state.particles[state.particles.tot - 1], state.particles);
        empty = {};
        e->update( state.particles.get_subset(state.particles.tot - 1), std::move(empty) );

        elDE *= constants::lB;
        nomN[scale] += sElDE * std::exp(-(sSrDE + elDE));
        denomN[scale] += std::exp(-(sSrDE + elDE));  
    
        for(int l = 0; l < state.particles.tot - 1; l++){
            state.particles.particles[l]->q = state.particles.particles[l]->q < 0.0 ? -1.0 : 1.0;
        }
    }
    nomN[0] += elDE * std::exp(-sSrDE);
    denomN[0] += std::exp(-sSrDE);

    //if(!state.overlap(state.particles.tot - 1)){
        facN += std::exp(-sSrDE);
    //}
    nSamples++;
    
    state.particles.remove(state.particles.tot - 1);

    this->samples++;
}

void Samplers::ModifiedWidomCoulomb::save(){
    double integralP = 0.0;
    for(int i = 0; i <= 10; i++){
        integralP += (nomP[i] / this->pSamples) / (denomP[i] / this->pSamples) * 1.0 / 11.0;
    }
    double integralN = 0.0;
    for(int i = 0; i <= 10; i++){
        integralN += (nomN[i] / this->nSamples) / (denomN[i] / this->nSamples) * 1.0 / 11.0;
    }


    double sumP = 0.5 * (nomP[0] / this->pSamples) / (denomP[0] / this->pSamples) + 0.5 * (nomP[10] / this->pSamples) / (denomP[10] / this->pSamples);
    double sumN = 0.5 * (nomN[0] / this->nSamples) / (denomN[0] / this->nSamples) + 0.5 * (nomN[10] / this->nSamples) / (denomN[10] / this->nSamples);

    printf("MW + sr: %lf el %lf tot %lf\n", -std::log(this->facP / this->pSamples), integralP, -std::log(this->facP / this->pSamples) + integralP);
    printf("MW - sr: %lf el %lf tot %lf\n", -std::log(this->facN / this->nSamples), integralN, -std::log(this->facN / this->nSamples) + integralN);
    printf("MWSum + el: %lf tot %lf\n", sumP, -std::log(this->facP / this->pSamples) + sumP);
    printf("MWSum - el: %lf tot %lf\n", sumN, -std::log(this->facN / this->nSamples) + sumN);
}

void Samplers::ModifiedWidomCoulomb::close(){};

/*
class Potential : public SamplerBase{
    private:

    double binWidth, dl;
    std::vector< std::vector<double> > potential;
    int bins;

    public:

    Potential(double dl, double binWidth){
        this->binWidth = binWidth;
        this->bins = dl / binWidth;
        this->potential.resize(this->bins, std::vector<double>(2));
        this->dl = dl;

        for(auto p : this->potential){
            p[0] = 0.0;
            p[1] = 0.0;
        }
    }

    void sample(Particles& particles){
        particles.add(T com, T pos, double r, 2.5, 1.0, 0.0, std::string name, bool image = false);
        this->potential[(int)particles.particles.back()->pos[2] + 0.5 * this->dl][0] += get_energy();
        this->potential[(int)particles.particles.back()->pos[2] + 0.5 * this->dl][1] += 1;
        particles.remove(particles.particles.back()->index);
    }

    void save(std::string filename){
        std::ofstream f ("potential_" + filename + ".txt");
        if (f.is_open())
        {
            for(unsigned int i = 0; i < this->potential.size(); i++){
                f << std::fixed << std::setprecision(10) << i * this->binWidth + this->binWidth / 2.0 << " " <<  
                     (double) this->potential[i][0] / this->potential[i][1] << "\n";
            }
            f.close();
        }
        else std::cout << "Unable to open file";
    }
};
*/

}