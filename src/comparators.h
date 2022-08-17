namespace Simulator{

namespace comparators{
    struct {
        template<typename T>
        bool operator()(T& a, T& b) const{  
            return a->weight < b->weight;
        }   
    } mLess;

    struct {
        bool operator()(const double a, Moves::Move* b) const{   
            return a < b->weight;
        }  

        bool operator()(Moves::Move* a, const double b) const{   
            return a->weight < b;
        }  
    } mdLess;
}

}