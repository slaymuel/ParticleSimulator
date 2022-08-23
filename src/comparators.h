namespace Simulator{

namespace comparators{
    struct {
        // Main template
        template<typename T>
        bool operator()(const T& a, const T& b) const{  
            return a->weight < b->weight;
        }
        template<typename T>
        bool operator()(const double a, const T& b) const{   
            return a < b->weight;
        }
        template<typename T>
        bool operator()(const T& a, const double b) const{   
            return a->weight < b;
        }  
    } mLess;
}

}