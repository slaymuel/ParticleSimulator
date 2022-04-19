#pragma once

//#include <source_location>
#include <chrono>
#include "logger.h"

namespace Simulator{

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
            std::chrono::duration<float> duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);

            //std::cout << location.function_name() << " took" << duration.count() << "s" << std::endl;
            Logger::Log(location, " took ", duration.count(), "s");
        }
    };

}