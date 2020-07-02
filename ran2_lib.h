#ifndef RAN_2_H
#define RAN_2_H

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

extern "C"{
    float ran2_(int*);
}

class ran2{
    public:
    static int ran_input;

    static double get_random(){
        float ran2_var = ran2_(&ran_input);
        //printf("rand: %lf\n", (double)ran2_var);
        return (double)ran2_var;
    }
};

int ran2::ran_input = -1*(int) time(NULL)*100*rand();
#endif
