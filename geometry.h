#include <vector>
#include <Eigen/Dense>

template <typename T>
class Geometry{
    public:

    std::vector<double> dh;  //half dimensions

    void resize(){
        static_cast<T*>(this)->resize_implementation();
    }

    bool is_inside(Eigen::VectorXd pos){
        static_cast<T*>(this)->is_inside_implementation(pos);
    }

};

class Cuboid : public Geometry<Cuboid>{
    void resize_implementation(){

    }

    bool is_inside_implementation(Eigen::VectorXd pos){
        if(pos[0] <= dh[0] && pos[0] >= -dh[0] &&
           pos[1] <= dh[1] && pos[1] >= -dh[1] &&
           pos[2] <= dh[2] && pos[2] >= -dh[2]){
               return true;
           }
        return false;
    }
};