#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "partialDifferential.h"

int main(int argc, char *argv[]) {
    int nx, ny, it, arg;
    char *outputFileName;
    FILE *outputFile = NULL;

    nx = ny = it = 0;

    for (arg = 1; arg < argc; arg++) {
        if (strcmp("-nx", argv[arg]) == 0) {
            arg++;
            nx = atoi(argv[arg]);
        }

        if (strcmp("-ny", argv[arg]) == 0) {
            arg++;
            ny = atoi(argv[arg]);
        }

        if (strcmp("-i", argv[arg]) == 0) {
            arg++;
            it = atoi(argv[arg]);
        }

        if (strcmp("-o", argv[arg]) == 0) {
            arg++;
            outputFileName = argv[arg];
            outputFile = fopen(outputFileName, "w");
        }
    }

    if (nx > 0 && ny > 0 && it > 0) {
        linearSystem linSys = initLinearSystem(nx, ny);

        setLinearSystem(&linSys);

        gaussSeidel(&linSys, it, outputFile);

        printOutput(&linSys, outputFile);

    } else {
        fprintf(stderr, "Argumentos incorretos. O formato deve ser: \"pdeSolver -nx <Nx> -ny <Ny> -i <maxIter> -o arquivo_saida\".\n");

        return -1;
    }

    return 0;
}
