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

        setLinearSystem(&linSys);

        gaussSeidel(&linSys, &it);

        // printf("Inferior afastada\n\n");
        // for (int i = 0; i < linSys.nx * linSys.ny; i++) {
        //     printf("%lf \n", linSys.iid[i]);
        // }

        // printf("Inferior\n\n");
        // for (int i = 0; i < linSys.nx * linSys.ny; i++) {
        //     printf("%lf \n", linSys.id[i]);
        // }

        // printf("Principal\n\n");
        // for (int i = 0; i < linSys.nx * linSys.ny; i++) {
        //     printf("%lf \n", linSys.md[i]);
        // }

        // printf("Superior\n\n");
        // for (int i = 0; i < linSys.nx * linSys.ny; i++) {
        //     printf("%lf \n", linSys.sd[i]);
        // }

        // printf("Superior superior\n\n");
        // for (int i = 0; i < linSys.nx * linSys.ny; i++) {
        //     printf("%lf \n", linSys.ssd[i]);
        // }

        // printf("Vetor b\n\n");
        // for (int i = 0; i < linSys.nx * linSys.ny; i++) {
        //     printf("%lf \n", linSys.b[i]);
        // }

        printf("Vetor solução.\n\n");
        for (int i = 0; i < linSys.nx * linSys.ny; i++) {
            printf("%lf\n", linSys.x[i]);
        }

    } else {
        fprintf(stderr, "Argumentos incorretos. O formato deve ser: \"pdeSolver -nx <Nx> -ny <Ny> -i <maxIter> -o arquivo_saida\".\n");

        return -1;
    }

    return 0;
}
