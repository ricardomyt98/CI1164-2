#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "partialDifferential.h"
#include "utils.h"

#define M_PI 3.14159265358979323846
#define SQR_PI M_PI *M_PI
#define X_Y_FUNCTION(i, j) (4 * SQR_PI) * ((sin(2 * M_PI * i) * sinh(M_PI * j)) + (sin(2 * M_PI * (M_PI - i)) * (sinh(M_PI * (M_PI - j)))))

/**
 * @brief Function to allocate space in memory.
 *
 * @param nx Number of points in x.
 * @param ny Number of points in y.
 * @return linearSystem Linear system struct.
 */
linearSystem initLinearSystem(int nx, int ny) {
    linearSystem linSys;

    linSys.ssd = malloc((nx * ny) * sizeof(real_t));
    linSys.sd = malloc((nx * ny) * sizeof(real_t));
    linSys.md = malloc((nx * ny) * sizeof(real_t));
    linSys.id = malloc((nx * ny) * sizeof(real_t));
    linSys.iid = malloc((nx * ny) * sizeof(real_t));

    linSys.b = malloc((nx * ny) * sizeof(real_t));

    linSys.x = malloc((nx * ny) * sizeof(real_t));

    linSys.nx = nx;
    linSys.ny = ny;

    return linSys;
}

/**
 * @brief Set the Linear System object
 *
 * @param linSys Linear system struct.
 */
void setLinearSystem(linearSystem *linSys) {
    real_t hx, hy, sqrHx, sqrHy;

    memset(linSys->x, 0.0, (linSys->nx * linSys->ny) * sizeof(real_t));  // Initial solution (zero).

    hx = M_PI / (linSys->nx + 1);
    hy = M_PI / (linSys->ny + 1);

    sqrHx = hx * hx;
    sqrHy = hy * hy;

    // ------------------------------------------------ FILL A DIAGONAL MATRIX ------------------------------------------------

    // Superior superior diagonal.
    for (int i = linSys->nx * linSys->ny; i > linSys->nx; i--) {
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

    for (int j = 1; j <= linSys->ny; j++) {
        for (int i = 1; i <= linSys->nx; i++) {
            linSys->b[idxB] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(0 + i * hx, 0 + j * hy);
            if (j == 1) {
                linSys->b[idxB] -= sin(2 * M_PI * (M_PI - i)) * sinh(SQR_PI);
            }

            if (i == 1) {
                linSys->b[idxB] -= 0;
            }

            if (j == linSys->ny) {
                linSys->b[idxB] -= sin(2 * M_PI * i) * sinh(SQR_PI);
            }

            if (i == linSys->nx) {
                linSys->b[idxB] -= 0;
            }

            idxB++;
        }
    }
}

/**
 * @brief
 *
 * @param linSys
 * @param res vetor residuo
 * @return real_t
 */
real_t l2Norm(linearSystem *linSys, real_t *res) {
    real_t *aux = malloc((linSys->nx * linSys->ny) * sizeof(real_t));

    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        aux[i] = linSys->b[i];
    }

    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        if (i - 1 >= 0) {
            aux[i] -= linSys->id[i] * linSys->x[i - 1];
        }

        if (i + 1 < linSys->nx * linSys->ny) {
            aux[i] -= linSys->sd[i] * linSys->x[i + 1];
        }

        if (i - linSys->nx >= 0) {
            aux[i] -= linSys->iid[i] * linSys->x[i - linSys->nx];
        }

        if (i + linSys->nx < linSys->nx * linSys->ny) {
            aux[i] -= linSys->ssd[i] * linSys->x[i + linSys->nx];
        }

        res[i] = aux[i];
    }

    real_t result = 0.0;

    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        result += aux[i] * aux[i];
    }

    return sqrt(result);
}

void printOutput(linearSystem *linSys, FILE *output) {
    int idxI, idxJ;
    real_t hx, hy;

    idxI = 0;
    idxJ = 0;
    hx = M_PI / (linSys->nx + 1);
    hy = M_PI / (linSys->ny + 1);

    fprintf(output, "---- Valores da função em cada ponto da grade ----\n");

    if (output) {
        for (int i = 0; i < linSys->nx * linSys->ny; i++) {
            fprintf(output, "(%lf, %lf) = %lf\n", idxI++ * hx, idxJ++ * hy, linSys->x[i]);
        }
    } else {
        for (int i = 0; i < linSys->nx * linSys->ny; i++) {
            printf("(%lf, %lf) = %lf\n", idxI++ * hx, idxJ++ * hy, linSys->x[i]);
        }
    }
}

void printGaussSeidelParameters(linearSystem *linSys, real_t *arrayItTime, real_t *arrayL2Norm, real_t *arrayResidue, FILE *output, int it) {
    real_t avrgTime = 0.0;

    fprintf(output, "###########\n");

    for (int i = 0; i < it; i++) {
        avrgTime += arrayItTime[i];
    }

    avrgTime /= it;

    fprintf(output, "# Tempo Método GS: %lfms\n", avrgTime);
    fprintf(output, "#\n");

    // ----------------------------------------------- Residue -----------------------------------------------
    fprintf(output, "# Resíduo\n");

    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        fprintf(output, "#i = %d : %lf\n", i, arrayResidue[i]);
    }

    fprintf(output, "#\n");

    // ----------------------------------------------- L2 Norm -----------------------------------------------

    fprintf(output, "# Norma L2 do Residuo\n");

    for (int i = 0; i < it; i++) {
        fprintf(output, "#i = %d : %lf\n", i, arrayL2Norm[i]);
    }

    fprintf(output, "###########\n\n");
}

/**
 * @brief Gauss Seidel function.
 *
 * @param linSys Linear system struct.
 * @param it Number of max iterations.
 */
void gaussSeidel(linearSystem *linSys, int *it, FILE *output) {
    real_t bk, itTime, *arrayL2Norm, *arrayItTime, *arrayResidue;

    int k = 1, i;
    arrayItTime = malloc(*it * sizeof(real_t));
    arrayL2Norm = malloc(*it * sizeof(real_t));
    arrayResidue = malloc((linSys->nx * linSys->ny) * sizeof(real_t));

    for (int k = 0; k < *(it); k++) {
        itTime = timestamp();
        for (int i = 0; i < linSys->nx * linSys->ny; i++) {
            bk = linSys->b[i];

            if (i - 1 >= 0) {
                bk -= linSys->id[i] * linSys->x[i - 1];
            }

            if (i + 1 < linSys->nx * linSys->ny) {
                bk -= linSys->sd[i] * linSys->x[i + 1];
            }

            if (i - linSys->nx >= 0) {
                bk -= linSys->iid[i] * linSys->x[i - linSys->nx];
            }

            if (i + linSys->nx < linSys->nx * linSys->ny) {
                bk -= linSys->ssd[i] * linSys->x[i + linSys->nx];
            }

            linSys->x[i] = bk / linSys->md[i];
        }

        arrayItTime[k] = timestamp() - itTime;
        arrayL2Norm[k] = l2Norm(linSys, arrayResidue);
    }

    printGaussSeidelParameters(linSys, arrayItTime, arrayL2Norm, arrayResidue, output, *it);
}