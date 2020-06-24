#pragma once

#include <vector>
#include <Eigen/Dense>
#include "particle.h"

class Geometry{
    private:

    public:

    std::vector<double> d;
    std::vector<double> _d;
    std::vector<double> dh;  //half dimensions
    std::vector<double> _dh;
    double volume;

    virtual void resize() = 0;
    virtual bool is_inside(std::shared_ptr<Particle>& p) = 0;
    virtual void pbc(std::shared_ptr<Particle>&& p) = 0;
    virtual double distance(Eigen::Vector3d& a, Eigen::Vector3d& b) = 0;
    virtual Eigen::Vector3d mirror(Eigen::Vector3d pos) = 0;
    virtual Eigen::Vector3d random_pos(double rf) = 0;
    virtual ~Geometry(){};

};


template<bool X = true, bool Y = true, bool Z = true>
class Cuboid : public Geometry{

    public:

    Cuboid(double x, double y, double z){
        this->d = {x, y, z};
        this->_d = this->d;
        this->dh = {x / 2.0, y / 2.0, z / 2.0};
        this->volume = x * y * z;
        printf("\tVolume: %lf\n", this->volume);
    }

    void resize(){

    }



    bool is_inside(std::shared_ptr<Particle>& p){
        //assert(p->com.size() == 3);

        if(X){
            if(p->com[0] >= this->dh[0] || p->com[0] <= -this->dh[0]){
                return false;
            }
        }
        else{
            if(p->com[0] + p->rf >= this->dh[0] || p->com[0] - p->rf <= -this->dh[0]){
                return false;
            }               
        }

        if(Y){
            if(p->com[1] >= this->dh[1] || p->com[1] <= -this->dh[1]){
                return false;
            }
        }
        else{
            if(p->com[1] + p->rf >= this->dh[1] || p->com[1] - p->rf <= -this->dh[1]){
                return false;
            }
        }

        if(Z){
            if(p->com[2] >= this->dh[2] || p->com[2] <= -this->dh[2]){
                return false;
            }
        }
        else{
            if(p->com[2] + p->rf >= this->dh[2] || p->com[2] - p->rf <= -this->dh[2]){
                return false;
            }
        }
        return true;
    }

    void pbc(std::shared_ptr<Particle>&& p){

        if(X){
            if(p->com[0] > this->dh[0]){
                p->com[0] -= d[0];
            }
            else if(p->com[0] < -this->dh[0]){
                p->com[0] += d[0];
            }
        }

        if(Y){
            if(p->com[1] > this->dh[1]){
                p->com[1] -= d[1];
            }
            else if(p->com[1] < -this->dh[1]){
                p->com[1] += d[1];
            }
        }

        if(Z){
            if(p->com[2] > this->dh[2]){
                p->com[2] -= d[2];
            }
            else if(p->com[2] < -this->dh[2]){
                p->com[2] += d[2];
            }
        }
        p->pos = p->com + p->qDisp;

        if(X){
            if(p->pos[0] > this->dh[0]){
                p->pos[0] -= d[0];
            }
            else if(p->pos[0] < -this->dh[0]){
                p->pos[0] += d[0];
            }
        }

        if(Y){
            if(p->pos[1] > this->dh[1]){
                p->pos[1] -= d[1];
            }
            else if(p->pos[1] < -this->dh[1]){
                p->pos[1] += d[1];
            }
        }

        if(Z){
            if(p->pos[2] > this->dh[2]){
                p->pos[2] -= d[2];
            }
            else if(p->pos[2] < -this->dh[2]){
                p->pos[2] += d[2];
            }
        }
    }

    double distance(Eigen::Vector3d& a, Eigen::Vector3d& b){
        Eigen::Vector3d disp = a - b;

        if(X){
            if(disp[0] > this->dh[0]){
                disp[0] -= this->d[0];
            }

            else if(disp[0] < -this->dh[0]){
                disp[0] += this->d[0];
            }
        }

        if(Y){
            if(disp[1] > this->dh[1]){
                disp[1] -= this->d[1];
            }

            else if(disp[1] < -this->dh[1]){
                disp[1] += this->d[1];
            }
        }

        if(Z){
            if(disp[2] > this->dh[2]){
                disp[2] -= this->d[2];
            }

            else if(disp[2] < -this->dh[2]){
                disp[2] += this->d[2];
            }
        }

        return disp.norm();
    }


    Eigen::Vector3d mirror(Eigen::Vector3d pos){
        Eigen::Vector3d m;
        m << pos[0], pos[1], (pos[2] >= 0) ? 2.0 * dh[2] - pos[2] : -2.0 * dh[2] - pos[2];
        return m;
    }

    Eigen::Vector3d random_pos(double rf){
        double x = (X) ? rf : 0.0;
        double y = (Y) ? rf : 0.0;
        double z = (Z) ? rf : 0.0;

        Eigen::Vector3d v;
        v = Random::get_vector();
        v << (dh[0] - x) * v[0], (dh[1] - y) * v[1], (dh[2] - z) * v[2];
        return v;
    }
};



template<bool X = true, bool Y = true, bool Z = true>
class CuboidImg : public Geometry{

    public:

    CuboidImg(double x, double y, double z){
        this->_d = {x, y,       z};
        this->d  = {x, y, 2.0 * z};
        this->dh = {this->d[0] / 2.0, this->d[1] / 2.0, this->d[2] / 2.0};
        this->_dh = {this->_d[0] / 2.0, this->_d[1] / 2.0, this->_d[2] / 2.0};

        this->volume = x * y * z;
        printf("\tBox dimensions: %.3lf, %.3lf, %.3lf\n", this->_d[0], this->_d[1], this->_d[2]);
        printf("\tVolume %.3lf\n", this->volume);
    }

    void resize(){}



    bool is_inside(std::shared_ptr<Particle>& p){
        if(p->com[0] >= this->_dh[0] || p->com[0] <= -this->_dh[0] ||
           p->com[1] >= this->_dh[1] || p->com[1] <= -this->_dh[1] ||
           p->com[2] + p->rf >= this->_dh[2] || p->com[2] - p->rf <= -this->_dh[2]){
               return false;
        }
        return true;
    }

    void pbc(std::shared_ptr<Particle>&& p){
        if(p->com[0] > this->_dh[0]){
            p->com[0] -= _d[0];
        }
        else if(p->com[0] < -this->_dh[0]){
            p->com[0] += _d[0];
        }


        if(p->com[1] > this->_dh[1]){
            p->com[1] -= _d[1];
        }
        else if(p->com[1] < -this->_dh[1]){
            p->com[1] += _d[1];
        }

        p->pos = p->com + p->qDisp;

        if(p->pos[0] > this->_dh[0]){
            p->pos[0] -= _d[0];
        }
        else if(p->pos[0] < -this->_dh[0]){
            p->pos[0] += _d[0];
        }


        if(p->pos[1] > this->_dh[1]){
            p->pos[1] -= _d[1];
        }
        else if(p->pos[1] < -this->_dh[1]){
            p->pos[1] += _d[1];
        }
    }

    double distance(Eigen::Vector3d& a, Eigen::Vector3d& b){
        Eigen::Vector3d disp = a - b;

        if(X){
            if(disp[0] > this->dh[0]){
                disp[0] -= this->d[0];
            }

            else if(disp[0] < -this->dh[0]){
                disp[0] += this->d[0];
            }
        }

        if(Y){
            if(disp[1] > this->dh[1]){
                disp[1] -= this->d[1];
            }

            else if(disp[1] < -this->dh[1]){
                disp[1] += this->d[1];
            }
        }


        if(Z){
            if(disp[2] > this->dh[2]){
                disp[2] -= this->d[2];
            }

            else if(disp[2] < -this->dh[2]){
                disp[2] += this->d[2];
            }
        }

        return disp.norm();
    }


    Eigen::Vector3d mirror(Eigen::Vector3d pos){
        Eigen::Vector3d m;
        m << pos[0], pos[1], (pos[2] >= 0) ? 2.0 * dh[2] - pos[2] : -2.0 * _dh[2] - pos[2];
        return m;
    }

    Eigen::Vector3d random_pos(double rf){
        Eigen::Vector3d v;
        v = Random::get_vector();
        v << _dh[0] * v[0], _dh[1] * v[1], (_dh[2] - rf) * v[2];
        return v;
    }
};



class Sphere : public Geometry{
    void resize(){

    }

    bool is_inside(std::shared_ptr<Particle>& p){
        if(sqrt(p->com[0] * p->com[0] + p->com[1] * p->com[1] + p->com[2] * p->com[2]) < this->dh[0]){
               return true;
           }
        return false;
    }

    double distance(Eigen::Vector3d& a, Eigen::Vector3d& b){ return 0.0; }

    void pbc(std::shared_ptr<Particle>&& p){  }

    Eigen::Vector3d mirror(Eigen::Vector3d pos){
        Eigen::Vector3d m;
        return m;
    }

    Eigen::Vector3d random_pos(double rf){
        Eigen::Vector3d v;
        return v;
    }
};