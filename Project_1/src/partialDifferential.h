#ifndef __PARTIAL_DIFFERENTIAL__
#define __PARTIAL_DIFFERENTIAL__

typedef double real_t;

typedef struct linearSystem {
    real_t *ssd;  // Superior superior diagonal.
    real_t *sd;   // Superior diagonal.
    real_t *md;   // Main diagonal.
    real_t *id;   // Inferior diagonal.
    real_t *iid;  // Inferior inferior diagonal.
    real_t *b;    // Independent terms.
    real_t *x;    // Solution.
    int nx, ny;
} linearSystem;

linearSystem initLinearSystem(int nx, int ny);

void setLinearSystem(linearSystem *linSys);

void gaussSeidel(linearSystem *linSys, int *it, FILE *output);

void printGaussSeidelParameters(real_t avrgTime, real_t *arrayL2Norm, FILE *output, int it);

real_t l2Norm(linearSystem *linSys);

void printOutput(linearSystem *linSys, FILE *output);

#endif  // __PARTIAL_DIFFERENTIAL__
