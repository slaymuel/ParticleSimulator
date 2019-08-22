#include <vector>
#include <Eigen/Dense>

class Geometry{
    private:

    public:

    std::vector<double> d;
    std::vector<double> dh;  //half dimensions

    virtual void resize() = 0;
    virtual bool is_inside(Eigen::Vector3d& pos) = 0;
    virtual double distance(Eigen::Vector3d& a, Eigen::Vector3d& b) = 0;
    virtual ~Geometry(){};

};


template<bool X, bool Y, bool Z>
class Cuboid : public Geometry{
    public:

    Cuboid(double x, double y, double z){
        this->d = {x, y, z};
        this->dh = {x / 2.0, y / 2.0, z / 2.0};
    }

    void resize(){

    }



    bool is_inside(Eigen::Vector3d& pos){
        assert(pos.size() == 3 && "Position is malformed...");

        if(X){
            if(pos[0] > this->dh[0]){
                pos[0] -= d[0];
            }
            else if(pos[0] < -this->dh[0]){
                pos[0] += d[0];
            }
        }
        if(Y){
            if(pos[1] > this->dh[1]){
                pos[1] -= d[1];
            }
            else if(pos[1] < -this->dh[1]){
                pos[1] += d[1];
            }
        }

        if(Z){
            if(pos[2] > this->dh[2]){
                pos[2] -= d[2];
            }
            else if(pos[2] < -this->dh[2]){
                pos[2] += d[2];
            }
        }

        if(pos[0] >= this->dh[0] || pos[0] <= -this->dh[0] ||
           pos[1] >= this->dh[1] || pos[1] <= -this->dh[1] ||
           pos[2] >= this->dh[2] || pos[2] <= -this->dh[2]){
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
};



class Sphere : public Geometry{
    void resize(){

    }

    bool is_inside(Eigen::Vector3d& pos){
        if(sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]) < this->dh[0]){
               return true;
           }
        return false;
    }

    double distance(Eigen::Vector3d& a, Eigen::Vector3d& b){ return 0.0; }
};