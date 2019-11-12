#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <likwid.h>

#include "partialDifferential.h"
#include "utils.h"

#define M_PI 3.14159265358979323846
#define SQR_PI M_PI *M_PI
#define X_Y_FUNCTION(i, j) (4 * SQR_PI) * ((sin(2 * M_PI * (i)) * sinh(M_PI * (j))) + (sin(2 * M_PI * (M_PI - (i))) * (sinh(M_PI * (M_PI - (j))))))

/**
 * @brief Function to allocate space in memory.
 *
 * @param nx Number of points in x.
 * @param ny Number of points in y.
 * @return linearSystem Linear system struct.
 */
linearSystem initLinearSystem(int nx, int ny) {
    linearSystem linSys;

    linSys.ssd = (real_t *)malloc((nx * ny) * sizeof(real_t));
    memset(linSys.ssd, 0.0, (nx * ny) * sizeof(real_t));

    linSys.sd = (real_t *)malloc((nx * ny) * sizeof(real_t));
    memset(linSys.sd, 0.0, (nx * ny) * sizeof(real_t));

    linSys.md = (real_t *)malloc((nx * ny) * sizeof(real_t));
    memset(linSys.md, 0.0, (nx * ny) * sizeof(real_t));

    linSys.id = (real_t *)malloc((nx * ny) * sizeof(real_t));
    memset(linSys.id, 0.0, (nx * ny) * sizeof(real_t));

    linSys.iid = (real_t *)malloc((nx * ny) * sizeof(real_t));
    memset(linSys.iid, 0.0, (nx * ny) * sizeof(real_t));

    linSys.b = (real_t *)malloc((nx * ny) * sizeof(real_t));

    linSys.x = (real_t *)malloc((nx * ny) * sizeof(real_t));
    memset(linSys.x, 0.0, (nx * ny) * sizeof(real_t));

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

    hx = M_PI / (linSys->nx + 1);
    hy = M_PI / (linSys->ny + 1);

    sqrHx = hx * hx;
    sqrHy = hy * hy;

    // ------------------------------------------------ FILL A DIAGONAL MATRIX ------------------------------------------------

    // Superior superior diagonal.
    for (int i = 0; i < (linSys->nx * linSys->ny) - linSys->nx; i++) {
        linSys->ssd[i] = sqrHx * (hy - 2);
    }

    // Superior diagonal.
    for (int k = 0; k < linSys->nx * linSys->ny; k++) {
        if ((k + 1) % linSys->nx != 0) {
            linSys->sd[k] = sqrHy * (hx - 2);
        }
    }

    // Main diagonal
    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        linSys->md[i] = 4 * (sqrHy + sqrHx + 2 * SQR_PI * sqrHx * sqrHy);
    }

    // Inferior diagonal
    for (int k = 0; k < linSys->nx * linSys->ny; k++) {
        if (k % linSys->nx != 0) {
            linSys->id[k] = sqrHy * (-2 - hx);
        }
    }

    // Inferior inferior diagonal
    for (int i = linSys->nx; i < linSys->nx * linSys->ny; i++) {
        linSys->iid[i] = sqrHx * (-2 - hy);
    }

    // ------------------------------------------------ FILL B ARRAY ------------------------------------------------

    int idxB = 0;

    for (int j = 1; j <= linSys->ny; j++) {
        for (int i = 1; i <= linSys->nx; i++) {
            linSys->b[idxB] = (2 * sqrHx * sqrHy) * X_Y_FUNCTION(0 + i * hx, 0 + j * hy);

            // Bottom edge of the mesh.
            if (j == 1) {
                linSys->b[idxB] -= (sin(2 * M_PI * (M_PI - i * hx)) * sinh(SQR_PI)) * (sqrHx * (-2 - hy));
            }

            // Left edge of the mesh (zero).
            if (i == 1) {
                linSys->b[idxB] -= 0;
            }

            // Upper edge of the mesh.
            if (j == linSys->ny) {
                linSys->b[idxB] -= (sin(2 * M_PI * i * hx) * sinh(SQR_PI)) * (sqrHx * (hy - 2));
            }

            // Right edge of the mesh (zero).
            if (i == linSys->nx) {
                linSys->b[idxB] -= 0;
            }

            idxB++;
        }
    }
}

/**
 * @brief Function to calculate L2 norm.
 *
 * @param linSys Linear system struct.
 * @return real_t
 */
real_t l2Norm(linearSystem *linSys) {
    real_t *aux = (real_t *)malloc((linSys->nx * linSys->ny) * sizeof(real_t));

    // Copying array "b" to an "aux" array.
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

        aux[i] -= linSys->md[i] * linSys->x[i];

        if (i - linSys->nx >= 0) {
            aux[i] -= linSys->iid[i] * linSys->x[i - linSys->nx];
        }

        if (i + linSys->nx < linSys->nx * linSys->ny) {
            aux[i] -= linSys->ssd[i] * linSys->x[i + linSys->nx];
        }
    }

    real_t result = 0.0;

    for (int i = 0; i < linSys->nx * linSys->ny; i++) {
        result += aux[i] * aux[i];
    }

    return sqrt(result);
}

/**
 * @brief Function to print the discretization matrix.
 *
 * @param linSys LinearSystem structure
 * @param output Output file.
 */
void printMesh(linearSystem *linSys, FILE *output) {
    real_t hx, hy;

    hx = M_PI / (linSys->nx + 1);
    hy = M_PI / (linSys->ny + 1);

    if (!output) {
        output = stdout;
    }

    int k = 0;

    for (int j = 1; j < linSys->ny; j++) {
        for (int i = 1; i < linSys->nx; i++) {
            fprintf(output, "%lf %lf %lf\n", i * hx, j * hy, linSys->x[k++]);
        }
    }
}

/**
 * @brief Function to print Gauss Seidel parameters.
 *
 * @param avrgTime Average time.
 * @param arrayL2Norm Array with all L2 norms.
 * @param output Output file.
 * @param it Number of iterations.
 */
void printGaussSeidelParameters(real_t avrgTime, real_t *arrayL2Norm, FILE *output, int it) {
    fprintf(output, "###########\n");

    fprintf(output, "# Tempo Método GS: %lfms\n", avrgTime);
    fprintf(output, "#\n");

    // ----------------------------------------------- L2 Norm -----------------------------------------------

    fprintf(output, "# Norma L2 do Residuo\n");

    for (int i = 0; i < it; i++) {
        fprintf(output, "#i = %d : %lf\n", i, arrayL2Norm[i]);
    }

    fprintf(output, "###########\n");
}

/**
 * @brief Gauss Seidel function.
 *
 * @param linSys Linear system struct.
 * @param it Number of max iterations.
 */
void gaussSeidel(linearSystem *linSys, int it, FILE *output) {
    real_t bk, itTime, *arrayL2Norm, acumItTime;

    acumItTime = 0.0;
    arrayL2Norm = (real_t *)malloc(it * sizeof(real_t));

    LIKWID_MARKER_START("Gauss_Seidel_Likwid_Performance");
    for (int k = 0; k < it; k++) {
        itTime = timestamp();
        for (int i = 0; i < linSys->nx * linSys->ny; i++) {
            bk = linSys->b[i];

            if (i - 1 >= 0) {
                bk -= linSys->id[i] * linSys->x[i - 1];
            }

            if (i + 1 < (linSys->nx * linSys->ny) - 1) {
                bk -= linSys->sd[i] * linSys->x[i + 1];
            }

            if (i - linSys->nx >= 0) {
                bk -= linSys->iid[i] * linSys->x[i - linSys->nx];
            }

            if (i + linSys->nx < (linSys->nx * linSys->ny) - 1) {
                bk -= linSys->ssd[i] * linSys->x[i + linSys->nx];
            }

            linSys->x[i] = bk / linSys->md[i];
        }

        acumItTime += timestamp() - itTime;
        LIKWID_MARKER_START("l2_Norm_Likwid_Performance");
        arrayL2Norm[k] = l2Norm(linSys);
        LIKWID_MARKER_STOP("l2_Norm_Likwid_Performance");
    }
    LIKWID_MARKER_STOP("Gauss_Seidel_Likwid_Performance");

    printGaussSeidelParameters(acumItTime / (it), arrayL2Norm, output, it);
}
