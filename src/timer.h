#pragma once

//#include <source_location>
#include <chrono>

struct Timer{
    private:
    std::chrono::steady_clock::time_point start;
    //std::source_location location;
    std::string location;

    public:
    Timer(std::string location){
        start = std::chrono::high_resolution_clock::now();
        this->location = location;
    }

    ~Timer(){
        std::chrono::duration<float> duration = std::chrono::high_resolution_clock::now() - start;

        //std::cout << location.function_name() << " took" << duration.count() << "s" << std::endl;
        std::cout << location<< " took " << duration.count() << "s" << std::endl;
    }
};