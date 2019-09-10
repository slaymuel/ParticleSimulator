namespace math{

    template<typename T, typename G>
    inline auto dot(T vec1, G vec2) -> double{
        return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    }


    template<typename T>
    inline T erfc_x( T x )
    {
        //static_assert(std::is_floating_point<T>::value, "type must be floating point");
        if(x < 0){
            return ( 2.0 - erfc_x(-x) );
        }
        T t = 1.0 / (1.0 + 0.3275911 * x);
        const T a1 = 0.254829592;
        const T a2 = -0.284496736;
        const T a3 = 1.421413741;
        const T a4 = -1.453152027;
        const T a5 = 1.061405429;
        return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) *  exp(-x * x);
        //return 1;
    }




    template<typename T>
    inline T erf_x( T x ) {
        return (1 - erfc_x(x));
    }





    template<typename T>
    inline auto norm(T x) -> double{
        double norm = 0;

        norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return std::sqrt(norm);
    }
}