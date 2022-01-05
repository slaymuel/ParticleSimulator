#include <cstdlib>
#include "../src/random/random.h"

int main(){
  double sum = 0;
  
  for(int i = 0; i < 10000000; i++){
    sum += 2.0*Simulator::Random::get_random() - 1.0;
  }
  
  std::cout << "Sum is " << sum << std::endl;
  
  if(std::abs(sum) < 2000.0) //Test Success
    return 0;
  else          //Test Failed
    return 1;
}
