#!/bin/sh

time gcc -c main.cpp particle.cpp -std=c++17
gcc -o main main.o particle.o -lstdc++