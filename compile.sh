gcc -c mormon.cpp -O3 -march=native -ffast-math -Xpreprocessor -fopenmp -std=c++17 -pthread -Wall
gcc -o mormon mormon.o -O3 -march=native -ffast-math -lstdc++ -lomp