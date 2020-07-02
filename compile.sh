#!/bin/sh

#gcc -c mormon.cpp Faddeeva.cpp xdrfile.c xdrfile_xtc.c xdrfile_trr.c -O3 -march=native -ffast-math -D DEBUG -Wall -Xpreprocessor -fopenmp -I/usr/local/include -L/usr/local/lib -lomp -std=c++17
#gcc -o mormon mormon.o Faddeeva.o xdrfile.o xdrfile_xtc.o xdrfile_trr.o -lstdc++ -O3 -ffast-math -Xpreprocessor -fopenmp -I/usr/local/include -L/usr/local/lib -lomp

#gfortran -c ran2.f
#gcc -c xdrfile.c xdrfile_xtc.c xdrfile_trr.c -O3 -march=native -ffast-math
#gcc -c mormon.cpp Faddeeva.cpp -O3 -march=native -ffast-math -D DEBUG -Wall -Xpreprocessor -fopenmp -I/usr/local/include -L/usr/local/lib -lomp -std=c++17
#gcc -o mormon mormon.o Faddeeva.o ran2.o xdrfile.o xdrfile_xtc.o xdrfile_trr.o -lstdc++ -O3 -ffast-math -Xpreprocessor -fopenmp -I/usr/local/include -L/usr/local/lib -lomp

#/usr/local/Cellar/gcc/9.3.0_1/bin/gfortran -c ran2.f
#/usr/local/Cellar/gcc/9.3.0_1/bin/gcc-9 -c xdrfile.c xdrfile_xtc.c xdrfile_trr.c -O3 -march=native -ffast-math
#/usr/local/Cellar/gcc/9.3.0_1/bin/gcc-9 -c mormon.cpp Faddeeva.cpp -O3 -march=native -ffast-math -D DEBUG -I/usr/local/include -fopenmp -std=c++17 -Wall -Wno-int-in-bool-context -Wno-reorder -Wno-sign-compare
#/usr/local/Cellar/gcc/9.3.0_1/bin/gcc-9 -o mormon mormon.o Faddeeva.o ran2.o xdrfile.o xdrfile_xtc.o xdrfile_trr.o -lstdc++ -O3 -ffast-math -lgomp -fopenmp

gfortran -c ran2.f
gcc -c xdrfile.c xdrfile_xtc.c xdrfile_trr.c -O2
gcc -c mormon.cpp Faddeeva.cpp -O2 -D DEBUG -Wall -Xpreprocessor -fopenmp -I/usr/local/include -L/usr/local/lib -lomp -std=c++17
gcc -o mormon mormon.o Faddeeva.o ran2.o xdrfile.o xdrfile_xtc.o xdrfile_trr.o -lstdc++ -O3 -ffast-math -Xpreprocessor -fopenmp -I/usr/local/include -L/usr/local/lib -lomp
