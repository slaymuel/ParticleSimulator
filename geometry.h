#include <vector>
#include <Eigen/Dense>

class Geometry{
    public:

    std::vector<double> dh;  //half dimensions

    virtual void resize() = 0;
    virtual bool is_inside(Eigen::VectorXd pos) = 0;
    virtual ~Geometry(){};

};

class Cuboid : public Geometry{
    public:
    Cuboid(double x, double y, double z){
        dh = {x, y, z};
    }
    void resize(){

    }

    bool is_inside(Eigen::VectorXd pos){
        assert(pos.size() == 3 && "Position is malformed...");
        if(pos[0] <= dh[0] && pos[0] >= -dh[0] &&
           pos[1] <= dh[1] && pos[1] >= -dh[1] &&
           pos[2] <= dh[2] && pos[2] >= -dh[2]){
               return true;
        }
        return false;
    }
};

class Sphere : public Geometry{
    void resize(){

    }

    bool is_inside(Eigen::VectorXd pos){
        if(sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]) < dh[0]){
               return true;
           }
        return false;
    }
};