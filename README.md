# GEM-pop
Implementing EM based mixture models for population genetic data. Three alternative algorithms is presented which can better circumvent the problem of multimodality (several solutions for the assignments of individuals into latent clusters).

## code compilation

compile command: 
gcc -o GEM GEM.c -Wall -I/usr/include -lm -lgsl -lgslcblas
