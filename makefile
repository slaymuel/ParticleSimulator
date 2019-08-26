CC = g++
CFLAGS = -O3 -Wall -std=c++14
all: 
	$(CC) $(CFLAGS) mormon.cpp -o mormon