#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "partialDifferential.h"

int main(int argc, char* argv[]) {
    double nx, ny;
    char* outputFile = NULL;
    int it, arg;

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
            outputFile = argv[arg];
        }
    }

    if (nx && ny && it && outputFile) {
        linearSystem linSys = initLinearSystem(nx, ny);

    } else {
        fprintf(stderr, "Argumentos incorretos. O formato deve ser: \"pdeSolver -nx <Nx> -ny <Ny> -i <maxIter> -o arquivo_saida\".\n");

        return -1;
    }

    return 0;
}
