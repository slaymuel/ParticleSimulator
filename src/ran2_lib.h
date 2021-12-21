#pragma once

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

extern "C"{
    float ran2_(int*);
}

class ran2{
    public:
    static int ran_input;
    static double get_random();
};