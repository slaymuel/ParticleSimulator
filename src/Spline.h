#include <vector>

class SplineData{

    public:
    std::vector<double> aKnots;
    std::vector<double> bKnots;
    std::vector<double >controlPoints;
    
    void load(std::vector<double> a, std::vector<double> b, std::vector<double >c){
        printf("\nLoading spline data.\n");
        this->aKnots = a;
        this->bKnots = b;
        this->controlPoints = c;
    }
};