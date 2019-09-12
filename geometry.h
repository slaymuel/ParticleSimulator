#pragma once

#include <vector>
#include <Eigen/Dense>
#include "particle.h"

class Geometry{
    private:

    public:

    std::vector<double> d;
    std::vector<double> dh;  //half dimensions
    double volume;

    virtual void resize() = 0;
    virtual bool is_inside(std::shared_ptr<Particle>& p) = 0;
    virtual double distance(Eigen::Vector3d& a, Eigen::Vector3d& b) = 0;
    virtual Eigen::Vector3d mirror(Eigen::Vector3d pos) = 0;
    virtual Eigen::Vector3d random_pos() = 0;
    virtual ~Geometry(){};

};


template<bool X = true, bool Y = true, bool Z = true>
class Cuboid : public Geometry{

    public:

    Cuboid(double x, double y, double z){
        this->d = {x, y, z};
        this->dh = {x / 2.0, y / 2.0, z / 2.0};
        this->volume = x * y * z;
    }

    void resize(){

    }



    bool is_inside(std::shared_ptr<Particle>& p){
        assert(p->pos.size() == 3);

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
        else{
            if(p->pos[2] + p->r > this->dh[2] || p->pos[2] - p->r < -this->dh[2]){
                return false;
            }
        }

        if(p->pos[0] >= this->dh[0] || p->pos[0] <= -this->dh[0] ||
           p->pos[1] >= this->dh[1] || p->pos[1] <= -this->dh[1] ||
           p->pos[2] >= this->dh[2] || p->pos[2] <= -this->dh[2]){
               return false;
        }
        return true;
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

    Eigen::Vector3d random_pos(){
        Eigen::Vector3d v;
        v = Random::get_vector();
        v << (dh[0] - 2.5) * (v[0] * 2.0 - 1), (dh[1] - 2.5) * (v[1] * 2.0 - 1), (dh[2] - 2.5) * (v[2] * 2.0 - 1);
        return v;
    }
};




class CuboidImg : public Geometry{
    private:
    std::vector<double> _dh;
    std::vector<double> _d;

    public:

    CuboidImg(double x, double y, double z){
        this->d = {x, y, 2.0 * z};
        this->dh = {this->d[0] / 2.0, this->d[1] / 2.0, this->d[2] / 2.0};
        this->_dh = {this->dh[0], this->dh[1], this->dh[2] / 2.0};
        this->_d  = {this->d[0],  this->d[1],  this->d[2] / 2.0};
        this->volume = x * y * z;
    }

    void resize(){

    }



    bool is_inside(std::shared_ptr<Particle>& p){
        assert(p->pos.size() == 3);

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


        if(p->pos[0] >= this->_dh[0] || p->pos[0] <= -this->_dh[0] ||
           p->pos[1] >= this->_dh[1] || p->pos[1] <= -this->_dh[1] ||
           p->pos[2] >= this->_dh[2] || p->pos[2] <= -this->_dh[2]){
               return false;
        }
        return true;
    }


    double distance(Eigen::Vector3d& a, Eigen::Vector3d& b){
        Eigen::Vector3d disp = a - b;


        if(disp[0] > this->dh[0]){
            disp[0] -= this->d[0];
        }

        else if(disp[0] < -this->dh[0]){
            disp[0] += this->d[0];
        }



        if(disp[1] > this->dh[1]){
            disp[1] -= this->d[1];
        }

        else if(disp[1] < -this->dh[1]){
            disp[1] += this->d[1];
        }



        if(disp[2] > this->dh[2]){
            disp[2] -= this->d[2];
        }

        else if(disp[2] < -this->dh[2]){
            disp[2] += this->d[2];
        }


        return disp.norm();
    }


    Eigen::Vector3d mirror(Eigen::Vector3d pos){
        Eigen::Vector3d m;
        m << pos[0], pos[1], (pos[2] >= 0) ? 2.0 * dh[2] - pos[2] : -2.0 * _dh[2] - pos[2];
        return m;
    }

    Eigen::Vector3d random_pos(){
        Eigen::Vector3d v;
        v = Random::get_vector();
        v << (dh[0] - 2.5) * (v[0] * 2.0 - 1), (dh[1] - 2.5) * (v[1] * 2.0 - 1), (_dh[2] - 2.5) * (v[2] * 2.0 - 1);
        return v;
    }
};



class Sphere : public Geometry{
    void resize(){

    }

    bool is_inside(std::shared_ptr<Particle>& p){
        if(sqrt(p->pos[0] * p->pos[0] + p->pos[1] * p->pos[1] + p->pos[2] * p->pos[2]) < this->dh[0]){
               return true;
           }
        return false;
    }

    double distance(Eigen::Vector3d& a, Eigen::Vector3d& b){ return 0.0; }

    Eigen::Vector3d mirror(Eigen::Vector3d pos){
        Eigen::Vector3d m;
        return m;
    }

    Eigen::Vector3d random_pos(){
        Eigen::Vector3d v;
        return v;
    }
};