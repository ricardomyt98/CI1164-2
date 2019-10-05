#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "partialDifferential.h"

#define M_PI 3.14159265358979323846
#define SQR_PI M_PI *M_PI

// #define I_JM_FUNCTION (hx, hy, sqrHx, sqrHy) - 1 * (sqrHx * (2 + hy))                      // [i, j - 1]
// #define IM_J_FUNCTION (hx, hy, sqrHx, sqrHy)(sqrHy * (hx - 2))                             // [i - 1, j]
// #define I_J_FUNCTION (hx, hy, sqrHx, sqrHy)(4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)  // [i, j]
// #define IP_J_FUNCTION (hx, hy, sqrHx, sqrHy) - 1 * (sqrHy * (2 + hx))                      // [i + 1, j]
// #define I_JP_FUNCTION (hx, hy, sqrHx, sqrHy) sqrHx *(hy - 2)                               // [i, j + 1]

#define X_Y_FUNCTION(i, j) (4 * SQR_PI) * ((sin(2 * M_PI * i) * sinh(M_PI * j)) + (sin(2 * M_PI * (M_PI - i)) * (sinh(M_PI * (M_PI - j)))))

linearSystem initLinearSystem(int nx, int ny) {
    linearSystem linSys;

    linSys.ssd = malloc((nx * ny) * sizeof(real_t));
    linSys.sd = malloc((nx * ny) * sizeof(real_t));
    linSys.md = malloc((nx * ny) * sizeof(real_t));
    linSys.id = malloc((nx * ny) * sizeof(real_t));
    linSys.iid = malloc((nx * ny) * sizeof(real_t));

    linSys.b = malloc(n(nx * ny) * sizeof(real_t));

    linSys.nx = nx;
    linSys.ny = ny;

    return linSys;
}

void setLinearSystem(linearSystem *linSys) {
    real_t fxy, hx, hy, sqrHx, sqrHy;

    hx = M_PI / (linSys->nx + 1);
    hy = M_PI / (linSys->ny + 1);

    sqrHx = hx * hx;
    sqrHy = hy * hy;

    // ------------------------------------------------ FILL A DIAGONAL MATRIX ------------------------------------------------

    // Superior superior diagonal.
    for (int i = linSys->nx * linSys->nx; i > linSys->nx; i--) {
        linSys->ssd[i] = 0;
    }
    for (int i = 0; i < (linSys->nx * linSys->ny) - linSys->nx; i++) {
        linSys->ssd[i] = sqrHx * (hy - 2);
    }

    // Superior diagonal.
    linSys->sd[linSys->nx * linSys->ny] = 0;
    for (int i = 0; i < (linSys->nx * linSys->ny) - 1; i++) {
        linSys->sd[i] = -1 * (sqrHy * (2 + hx));
    }

    // Main diagonal
    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        linSys->md[i] = (4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx);
    }

    // Inferior diagonal
    linSys->id[0] = 0;
    for (int i = 1; i < linSys->nx * linSys->ny; i++) {
        linSys->id[i] = sqrHy * (hx - 2);
    }

    // Inferior inferior diagonal
    for (int i = 0; i < linSys->nx; i++) {
        linSys->iid[i] = 0;
    }
    for (int i = linSys->nx; i < linSys->nx * linSys->ny; i++) {
        linSys->iid[i] = -1 * (sqrHx * (2 + hy));
    }

    // ------------------------------------------------ FILL B ARRAY ------------------------------------------------

    int idxB = 0;

    for (int j = 1; j < linSys->ny; j++) {
        for (int i = 1; i < linSys->nx; i++) {
            if (j == 1 && (i != 1 && i != linSys->ny)) {
                linSys->b[idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(0 + i * hx, 0 + j * hy);
            }

            if (j == linSys->nx && (i != 1 && i != linSys->ny)) {
                linSys->b[idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(0 + i * hx, 0 + j * hy);
            }

            if (i == 1) {
                linSys->b[idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(0 + i * hx, 0 + j * hy)) - (sin(2 * M_PI * (M_PI - i)) * sinh(SQR_PI));  // Para a primeira linha, b sempre será subtraido uma mesms contante (U(x, 0)).
            }

            if (i == linSys->ny) {
                linSys->b[idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(0 + i * hx, 0 + j * hy)) - (sin(2 * M_PI * i) * senh(SQR_PI));  // Para a última linha, b sempre será subtraido uma mesms contante (U(x, M_PI)).
            }
        }
    }
}
