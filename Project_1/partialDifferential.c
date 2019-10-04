#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "partialDifferential.h"

#define M_PI 3.14159265358979323846
#define SQR_PI M_PI *M_PI

#define HIGHER_VALUE(x, y) x > y ? x : y
#define LESSER_VALUE(x, y) x < y ? x : y

#define X_Y_FUNCTION(lin, col) (4 * SQR_PI) * ((sin(2 * M_PI * lin) * sinh(M_PI * col)) + (sin(2 * M_PI * (M_PI - lin)) * (sinh(M_PI * (M_PI - col)))))

linearSystem initLinearSystem(int nx, int ny) {
    linearSystem linSys;

    linSys.ssd = malloc(LESSER_VALUE(nx, ny) * sizeof(real_t));
    linSys.sd = malloc(LESSER_VALUE(nx, ny) * sizeof(real_t));
    linSys.md = malloc(HIGHER_VALUE(nx, nx) * sizeof(real_t));
    linSys.id = malloc(LESSER_VALUE(nx, ny) * sizeof(real_t));
    linSys.iid = malloc(LESSER_VALUE(nx, ny) * sizeof(real_t));

    linSys.b = malloc(ny * sizeof(real_t));

    linSys.nx = nx;
    linSys.ny = ny;

    linSys.idxSsd = 0;
    linSys.idxSd = 0;
    linSys.idxMd = 0;
    linSys.idxId = 0;
    linSys.idxIid = 0;
    linSys.idxB = 0;

    return linSys;
}

void setLinearSystem(linearSystem *linSys) {
    real_t fxy, hx, hy, sqrHx, sqrHy;

    hx = M_PI / (linSys->nx + 1);
    hy = M_PI / (linSys->ny + 1);

    sqrHx = hx * hx;
    sqrHy = hy * hy;

    for (int lin = 1; lin < linSys->nx; lin++) {
        for (int col = 1; col < linSys->nx; col++) {
            if (lin + 2 == col) {  // Superior superior
                if (lin == 1 && col == 3) {
                    linSys->ssd[linSys->idxSsd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->ssd[linSys->idxSsd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->ssd[linSys->idxSsd] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->ssd[linSys->idxSsd++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col)) - (sin(2 * M_PI * (M_PI - lin)) * sinh(SQR_PI));
                } else if (lin == linSys->nx - 1 && col == linSys->ny - 3) {
                    linSys->ssd[linSys->idxSsd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->ssd[linSys->idxSsd] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->ssd[linSys->idxSsd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->ssd[linSys->idxSsd++] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                         // [i + 1, j]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);
                } else {                                                                                              // middle
                    linSys->ssd[linSys->idxSsd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->ssd[linSys->idxSsd] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->ssd[linSys->idxSsd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->ssd[linSys->idxSsd] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->ssd[linSys->idxSsd++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);
                }
            } else if (lin + 1 == col) {  // superior
                if (lin == 1 && col == 2) {
                    linSys->sd[linSys->idxSd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->sd[linSys->idxSd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->sd[linSys->idxSd] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->sd[linSys->idxSd++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col)) - (sin(2 * M_PI * (M_PI - lin)) * sinh(SQR_PI));

                } else if (lin == linSys->nx - 1 && col == linSys->ny - 2) {
                    linSys->sd[linSys->idxSd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->sd[linSys->idxSd] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->sd[linSys->idxSd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->sd[linSys->idxSd++] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                         // [i + 1, j]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);

                } else {                                                                                            // middle
                    linSys->sd[linSys->idxSd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->sd[linSys->idxSd] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->sd[linSys->idxSd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->sd[linSys->idxSd] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->sd[linSys->idxSd++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);
                }
            } else if (lin == col) {  // main
                if (lin == 1 && col == 1) {
                    linSys->md[linSys->idxMd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->md[linSys->idxMd] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->md[linSys->idxMd++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col)) - (sin(2 * M_PI * (M_PI - lin)) * sinh(SQR_PI));

                } else if (lin == linSys->nx - 1 && col == linSys->ny - 1) {
                    linSys->md[linSys->idxMd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->md[linSys->idxMd] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->md[linSys->idxMd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->md[linSys->idxMd++] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                         // [i + 1, j]

                    linSys->b[linSys->idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col)) - (sin(2 * M_PI * lin) * senh(SQR_PI));

                } else {                                                                                            // middle
                    linSys->md[linSys->idxMd] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->md[linSys->idxMd] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->md[linSys->idxMd] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->md[linSys->idxMd] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->md[linSys->idxMd++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);
                }
            } else if (lin == col + 1) {  // inferior
                if (lin == 2 && col == 1) {
                    linSys->id[linSys->idxId] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->id[linSys->idxId] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->id[linSys->idxId] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->id[linSys->idxId++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);

                } else if (lin == linSys->nx - 2 && col == linSys->ny - 1) {
                    linSys->id[linSys->idxId] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->id[linSys->idxId] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->id[linSys->idxId] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->id[linSys->idxId++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col)) - (sin(2 * M_PI * lin) * senh(SQR_PI));

                } else {                                                                                            // middle
                    linSys->id[linSys->idxId] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->id[linSys->idxId] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->id[linSys->idxId] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->id[linSys->idxId] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->id[linSys->idxId++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);
                }
            } else if (lin == col + 2) {  // inferior inferior
                if (lin == 3 && col == 1) {
                    linSys->iid[linSys->idxIid] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->iid[linSys->idxIid] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->iid[linSys->idxIid] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->iid[linSys->idxIid++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);

                } else if (lin == linSys->ny - 1 && col == linSys->nx - 3) {
                    linSys->iid[linSys->idxIid] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->iid[linSys->idxIid] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->iid[linSys->idxIid] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->iid[linSys->idxIid++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = ((2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col)) - (sin(2 * M_PI * lin) * senh(SQR_PI));

                } else {                                                                                              // middle
                    linSys->iid[linSys->idxIid] = (-sqrHx * (2 + hy)) * (U[lin][col - 1]);                            // [i, j - 1]
                    linSys->iid[linSys->idxIid] += (sqrHy * (hx - 2)) * (U[lin - 1][col]);                            // [i - 1, j]
                    linSys->iid[linSys->idxIid] += ((4 * sqrHy) * (1 + sqrHx + 2 * SQR_PI * sqrHx)) * (U[lin][col]);  // [i, j]
                    linSys->iid[linSys->idxIid] += (-sqrHy * (2 + hx)) * (U[lin + 1][col]);                           // [i + 1, j]
                    linSys->iid[linSys->idxIid++] += (sqrHx * (hy - 2)) * (U[lin][col + 1]);                          // [i, j + 1]

                    linSys->b[linSys->idxB++] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(lin, col);
                }
            }
        }
    }
}
