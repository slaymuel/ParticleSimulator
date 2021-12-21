#include "random.h"

void Random::create_CDF(std::string filename){
    
    
    std::cout << "Creating CDF" << std::endl;
    std::ifstream f(filename);
    std::string line;
    double integral = 0;
    
    while (std::getline(f, line)) {
        double x, y;
        
        std::istringstream ss(line);
        ss >> x >> y;
        std::cout << x << " " << y << std::endl;

        if(y > 1e-10){
            customVals.push_back( std::vector<double>({x, y}) );
            integral += y;
        }
    }
    double dx = customVals[1][0] - customVals[0][0];
    std::cout << "dx is: " << dx << std::endl;
    //integral *= dx;
    for(std::vector<double>& val : customVals){
        val[1] /= integral;
    }

    for(int i = 1; i < customVals.size(); i++){
        customVals[i][1] += customVals[i - 1][1];
        printf("vals %lf %lf\n", customVals[i][0], customVals[i][1]);
    }
}

double Random::get_random_from_distribution(){
    double r = get_random();

    for(int i = 0; i < customVals.size(); i++){
        if(r < customVals[i][1]){
            return customVals[i][0];
        }
    }
    return -1.0;
}

double Random::get_random(){
    //return (*real_dist)(rand_gen);
    return ran2::get_random();
}

int Random::get_random(int i){
    return i * get_random();
}

Eigen::Vector3d Random::random_pos_box(double rf, std::vector<double> box){
    Eigen::Vector3d v;
    v = get_vector();
    v << box[0] * v[0], box[1] * v[1], (box[2] - rf) * v[2];
    return v;
}

Eigen::Vector3d Random::get_random_vector(double l){
    Eigen::Vector3d v(get_random()*2.0*l - l, get_random()*2.0*l - l, get_random()*2.0*l - l);
    if(l < 1e-10){
        return v;
    }
    while(v.norm() > l){
        v << get_random()*2.0*l - l, get_random()*2.0*l - l, get_random()*2.0*l - l;
    }
    return v;
}

Random::get_vector::operator Eigen::Vector3d(){
    Eigen::Vector3d v(get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0);
    return v;
}
Random::get_vector::operator std::vector<double>(){
    std::vector<double> v = {get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0, get_random() * 2.0 - 1.0};
    return v;
}

Random::get_norm_vector::operator Eigen::Vector3d(){
    double phi = get_random() * 2.0 * PI;
    double z = get_random() * 2.0 - 1.0;
    Eigen::Vector3d v(std::sqrt(1.0 - z*z) * std::cos(phi), std::sqrt(1.0 - z*z) * std::sin(phi), z);
    
    return v;
}
Random::get_norm_vector::operator std::vector<double>(){
    double phi = get_random() * 2.0 * PI;
    double z = get_random() * 2.0 - 1.0;
    std::vector<double> v = {std::sqrt(1.0 - z*z) * std::cos(phi), std::sqrt(1.0 - z*z) * std::sin(phi), z};
    
    return v;
}
