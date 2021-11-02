#!/bin/bash
gcc -Wall -c snu.c -std=c99
g++ -Wall -c Probability.cc
g++ -Wall Probability.o snu.o -lglobes -lgsl -lgslcblas -O3 -o Probability
./Probability 
rm *.o
#gnuplot plotter.gnuplot
