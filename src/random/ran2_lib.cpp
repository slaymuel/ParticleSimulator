#include "ran2_lib.h"

namespace Simulator{

double ran2::get_random(){
    float ran2_var = ran2_(&ran_input);
    return (double)ran2_var;
}

int ran2::ran_input = -1*(int) time(NULL)*100*rand();

}