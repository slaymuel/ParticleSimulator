#pragma once

//#include <source_location>
#include <chrono>
#include "logger.h"

namespace Simulator{

    // RAII timer
    struct Timer{
        private:
        std::chrono::high_resolution_clock::time_point start;
        std::string location;

        public:
        Timer(std::string location){
            start = std::chrono::high_resolution_clock::now();
            this->location = location;
        }

        ~Timer(){
            std::chrono::duration<float> duration = 
                    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
            Logger::Log(location, " took ", duration.count(), "s");
        }

        // Time a function 
        template<typename T>
        static void time(T& func){
            Timer timer(__FUNCTION__);
            func();
        }
    };

    #ifdef _TIMERS_
        #define TIMEIT Timer timer(__FUNCTION__);
    #else
        #define TIMEIT;
    #endif
}