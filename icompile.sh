icc -c mormon.cpp -O3 -march=native -ffast-math -fopenmp -std=c++14 -Wall
icc -o mormon mormon.o -O3 -march=native -ffast-math -lstdc++ -lomp