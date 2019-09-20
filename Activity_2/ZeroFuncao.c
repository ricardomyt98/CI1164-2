#include <math.h>
#include <stdio.h>

#include "ZeroFuncao.h"
#include "utils.h"

#define ENSURE_CONVERGENCE_CODE -1
#define ENSURE_CONVERGENCE_MSGE "Error, function results need to be opposed to ensure convergence.\n"

#define ZERO_DIVISION_CODE -2
#define ZERO_DIVISION_MSGE "Error, division by zero.\n"

#define MAX_IT_CODE -3
#define MAX_IT_MSGE "Error, the maximum number of allowed iterations reached.\n"

#define FUNCTION_RESULT_EQUAL_CODE -4
#define FUNCTION_RESULT_EQUAL_MSGE "The values of the results of fx0 and fx1 should be different to ensure convergence.\n"

#define POLINOMIAL_DEGREE_CODE -5
#define POLINOMIAL_DEGREE_MSGE "the polynomial degree must be greater than or equal to one.\n"

/**
 * @brief Bisection method.
 *
 * @param f Pointer to the function.
 * @param a Inferior limit.
 * @param b Upper limit.
 * @param eps Maximum relative error at x.
 * @param it Number of iterations.
 * @param raiz Root found.
 * @return int Return code.
 */
int bisseccao(double (*f)(const double x), double a, double b, double eps, int* it, double* raiz) {
    if (f(a) * f(b) >= 0) {
        fprintf(stderr, ENSURE_CONVERGENCE_MSGE);
        return ENSURE_CONVERGENCE_CODE;  // (Error) Opposite sings (requirement).
    } else {
        double rPrev, fa = f(a), fb = f(b), fr, rootFabs, relError;

        *raiz = (a + b) / 2;
        *it = 1;

        do {
            rPrev = *raiz;
            fr = f(rPrev);

            if (fa * fr < 0.0) {
                b = *raiz;  // Between a and raiz.
                fb = fr;
            } else if (fa * fr > 0.0) {
                a = *raiz;  // Between b and raiz.
                fa = fr;
            }

            *raiz = (a + b) / 2;
            rootFabs = fabs(*raiz);  // Absolute value of root.

            if (rootFabs == 0.0) {
                fprintf(stderr, ZERO_DIVISION_MSGE);
                return ZERO_DIVISION_CODE;  // (Error) Previne the division by zero to get the approximate error.
            } else {
                relError = fabs((*raiz - rPrev) / (*raiz)) * 100;
            }

            (*it)++;  // Increment the number of iterations.

            if ((*it) > MAXIT) {
                fprintf(stderr, MAX_IT_MSGE);
                return MAX_IT_CODE;
            }
        } while (relError > eps);  // Is not greater than the max number of iterations and calculate the approximate error.

#ifdef __DEBUG__
        fprintf(stderr, "\n\n****** Iterações Bissecção *******\n");
        fprintf(stderr, "=================================================================================================\n");
        fprintf(stderr, "  it            a            b           xs        f(xs)        |a-b|         f(a)         f(b)\n");
        fprintf(stderr, "=================================================================================================\n");
#endif

#ifdef __DEBUG__
        fprintf(stderr, "%4d %12.8g %12.8g %12.8g ", *it, a, b, *raiz);
        fprintf(stderr, "%12.4e %12.4e %12.4e %12.4e\n", fr, fabs(a - b), fa, fb);
#endif

        return 0;
    }
}

/**
 * @brief Newton method.
 *
 * @param f Pointer to the function.
 * @param df Pointer to the derivative of a function.
 * @param x0 Initial value.
 * @param eps Maximum relative error at x.
 * @param it Number of iterations.
 * @param raiz Root found.
 * @return int Return code.
 */
int newton(double (*f)(const double x), double (*df)(const double x), double x0, double eps, int* it, double* raiz) {
    if (df(x0) == 0.0) {
        fprintf(stderr, ZERO_DIVISION_MSGE);
        return ZERO_DIVISION_CODE;  // (Error) Division by zero of the initial guess.
    } else {
        double rPrev, dfPrev;

        *raiz = x0 - (f(x0) / df(x0));
        *it = 0;

        while (fabs((*raiz - rPrev) / (*raiz)) * 100 > eps) {
            rPrev = *raiz;
            dfPrev = df(rPrev);

            if (dfPrev == 0.0) {
                fprintf(stderr, ZERO_DIVISION_MSGE);
                return ZERO_DIVISION_CODE;  // (Error) Division by zero of the previous root.
            } else {
                *raiz = rPrev - (f(rPrev) / dfPrev);
            }

            (*it)++;  // Increment the number of iterations.

            if ((*it) == MAXIT) {
                fprintf(stderr, MAX_IT_MSGE);
                return MAX_IT_CODE;  // (Error)The number of iterations is not enough.
            }
        }

#ifdef __DEBUG__
        fprintf(stderr, "\n\n****** Iterações Newton-Raphson *******\n");
        fprintf(stderr, "=================================================================================================\n");
        fprintf(stderr, "  it           x0           x0         raiz       f(raiz)   |raiz-x0|        f(x0)       df(x0)\n");
        fprintf(stderr, "=================================================================================================\n");
#endif

#ifdef __DEBUG__
        fprintf(stderr, "%4d %12.8g %12.8g %12.8g ", *it, x0, x0, *raiz);
        fprintf(stderr, "%12.4e %12.4e %12.4e %12.4e\n", f(*raiz), fabs(*raiz - x0), fx, dfx);
#endif

        return 0;
    }
}

/**
 * @brief Secant method.
 *
 * @param f Pointer to the function.
 * @param x0 Initial value.
 * @param x1 Initial value.
 * @param eps Maximum relative error at x.
 * @param it Number of iterations.
 * @param raiz Root found.
 * @return int Return code.
 */
int secante(double (*f)(const double x), double x0, double x1, double eps, int* it, double* raiz) {
    if (f(x0) == f(x1)) {
        fprintf(stderr, FUNCTION_RESULT_EQUAL_MSGE);
        return FUNCTION_RESULT_EQUAL_CODE;  // (Error) The resulting values of the intervals in f() must be different to find a solution.
    } else if (f(x0) * f(x1) >= 0.0) {
        fprintf(stderr, ENSURE_CONVERGENCE_MSGE);
        return ENSURE_CONVERGENCE_CODE;  // (Error) Opposite sings (requirement).
    } else {
        *it = 1;

        do {
            *raiz = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0));
            x0 = x1;
            x1 = *raiz;

            (*it)++;  // Increment the number of iterations.

            if ((*it) == MAXIT) {
                fprintf(stderr, MAX_IT_MSGE);
                return MAX_IT_CODE;
            }
        } while (fabs(f(*raiz)) > eps);
    }

#ifdef __DEBUG__
    fprintf(stderr, "\n\n****** Iterações Secante *******\n");
    fprintf(stderr, "=================================================================================================\n");
    fprintf(stderr, "  it           x0           x1         raiz      f(raiz)    |raiz-x1|        f(x0)        f(x1)\n");
    fprintf(stderr, "=================================================================================================\n");
#endif

#ifdef __DEBUG__
    fprintf(stderr, "%4d %12.8g %12.8g %12.8g ", *it, x0, x1, *raiz);
    fprintf(stderr, "%12.4e %12.4e %12.4e %12.4e\n", f(*raiz), fabs(*raiz - x0), fx0, fx1);
#endif

    return 0;
}

/**
 * @brief Derivative and polynomial calculation.
 *
 * @param pol Polynomial structure.
 * @param x Value to calculate the derivative and polynomial.
 * @param px Polynomial result.
 * @param dpx Result of the derivative.
 * @return int Return code.
 */
int calcPolinomioEDerivada(Polinomio pol, double x, double* px, double* dpx) {
    if (pol.grau < 1) {
        fprintf(stderr, POLINOMIAL_DEGREE_MSGE);
        return POLINOMIAL_DEGREE_CODE;
    } else {
        double b = pol.p[pol.grau], c = b;

        for (int i = pol.grau - 1; i; i--) {
            b = pol.p[i] + b * x;
            c = b + c * x;
        }

        b = pol.p[0] + b * x;
        *px = b;
        *dpx = c;

        return 0;
    }
}

/**
 * @brief Average calculation.
 *
 * Possíveis erros de precisão:
 * Devido a, possivelmente, haver uma grande variação entre os valores do vetor, a soma destes pode gerar propagação de erro a cada iteração.
 * Uma maneira de minimizar este erro seria, para cada valor do vetor, ir dividindo este por "n" ("padronizando" os valores para a soma), porém isto geraria o dobro de operações, podendo causar grande lentidão.
 *
 * Possível erro de divisão:
 * A soma dos elementos do vetor pode resultar em um valor muito grande, "n" ser muito pequeno e vice-versa, causando erro na divisão.
 *
 * @param valores Values array.
 * @param n Array length.
 * @return double Average calculation result.
 */
double media(double* valores, unsigned long n) {
    double soma = 0.0;

    for (int i = 0; i < n; i++) {
        soma = soma + valores[i];
    }

    return soma / n;
}
