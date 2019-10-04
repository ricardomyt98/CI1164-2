#ifndef __PARTIAL_DIFFERENTIAL__
#define __PARTIAL_DIFFERENTIAL__

typedef float real_t;

typedef struct linearSystem {
    real_t *ssd;  // Superior superior diagonal.
    real_t *sd;   // Superior diagonal.
    real_t *md;   // Main diagonal.
    real_t *id;   // Inferior diagonal.
    real_t *iid;  // Inferior inferior diagonal.
    real_t *b;    // Independent terms.
    int nx, ny, idxSsd, idxSd, idxMd, idxId, idxIid, idxB;
} linearSystem;

linearSystem initLinearSystem(int nx, int ny);

void setLinearSystem(linearSystem *linSys);

#endif  // __PARTIAL_DIFFERENTIAL__
