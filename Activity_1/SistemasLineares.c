#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

#define DIV_ZERO_CODE -1
#define DIV_ZERO_MSGE "Error, division by zero.\n"

#define MAX_IT_CODE -2
#define MAX_IT_MSGE "Error, the maximum number of allowed iterations reached.\n"

/**
 * @brief Esta função calcula a norma L2 do resíduo de um sistema linear.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param x Solução do sistema linear.
 * @return real_t Final result.
 */
real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  real_t *tmp = malloc(SL->n * sizeof(real_t));

  for (int i = 0; i < SL->n; i++)
  {

    tmp[i] = SL->b[i];

    for (int j = 0; j < SL->n; j++)
    {
      tmp[i] -= SL->A[i * SL->n + j] * x[j]; // b – A * x0.
    }
  }

  real_t result = 0;

  for (int i = 0; i < SL->n; i++)
  {
    result += tmp[i] * tmp[i];
  }

  free(tmp);

  return result;
}

/**
 * @brief Function to find the maximum value of a column and return the index.
 *
 * @param SL Linear system (struct).
 * @param col Column index.
 * @return int The index of the highest value in a column.
 */
int encontraMax(SistLinear_t *SL, int col)
{
  int iPivo = col;

  for (int i = col + 1; i < SL->n; i++)
  {
    if (SL->A[i * SL->n + (i + 1)] > iPivo)
    {
      iPivo = i;
    }
  }

  return iPivo;
}

/**
 * @brief Function to change two position lines.
 *
 * @param SL Linear system (struct).
 * @param iPivo Number referred to pivo index.
 * @param currCol Current index number.
 */
void trocaPosicao(SistLinear_t *SL, int iPivo, int currCol)
{
  int aux;

  for (int i = 0; i < SL->n; i++)
  {
    aux = SL->A[i];
    SL->A[i] = SL->A[i + SL->n];
    SL->A[i + SL->n] = aux;
    aux = SL->b[i];
    SL->b[i] = SL->b[i + SL->n];
    SL->b[i + SL->n] = aux;
  }
}

/**
 * @brief Retro-replacement function.
 *
 * @param SL Linear system (struct).
 * @param x Solution array.
 */
void resolucaoRetroativa(SistLinear_t *SL, real_t *x)
{
  for (int i = SL->n - 1; i >= 0; i--)
  {
    real_t tmp = SL->b[i];

    for (int j = i + 1; j < SL->n; j++)
    {
      tmp -= SL->A[i * SL->n + (i + 1)] * x[i];
    }

    x[i] = tmp / SL->A[i * SL->n + i];
  }
}

/**
 * @brief Método da Eliminação de Gauss.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param x Ponteiro para o vetor solução.
 * @param pivotamento Flag para indicar se o pivotamento parcial deve ser feito (!=0).
 * @return Código de erro. 0 em caso de sucesso.
 */
int eliminacaoGauss(SistLinear_t *SL, real_t *x, int pivotamento)
{
  for (int i = 0; i < SL->n; i++)
  {
    if (pivotamento != 0)
    {
      int iPivo = encontraMax(SL, i); // Find the highest value and return the index.

      if (SL->A[iPivo * SL->n + i] == 0) // Likely not satisfied condition (floating point).
      {
        fprintf(stderr, DIV_ZERO_MSGE);

        return DIV_ZERO_CODE;
      }

      if (i != iPivo)
      {
        trocaPosicao(SL, iPivo, i); // Change lines.
      }

      for (int j = i + 1; j < SL->n; j++)
      {
        real_t m = SL->A[j * SL->n + i] / SL->A[i * SL->n + i]; // Calculate the variable "m" to multiply by the next current line.

        SL->A[j * SL->n + i] = 0.0;

        for (int k = i + 1; k < SL->n; k++)
        {
          SL->A[j * SL->n + k] -= SL->A[i * SL->n + k] * m;
        }

        SL->b[j] -= SL->b[i] * m;
      }
    }
  }

  resolucaoRetroativa(SL, x);
  double normaL2 = normaL2Residuo(SL, x);

  return (0);
}

/**
 * @brief Function to calculate the maximum alpha used in the Gauss Jacobi method.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param line Current line.
 * @return int Alpha value.
 */
int encontraMaxAlpha(SistLinear_t *SL, int line)
{
  real_t lastAlpha, alpha = 0;
  int jAlpha = 0; // Index of the column with the highest alpha.

  for (int j = line; j < SL->n; j++) // Column changes.
  {
    lastAlpha = alpha;
    alpha = 0;

    for (int k = line; k < j; k++) // Before the current target (matrix number).
    {
      alpha += fabs(SL->A[line * SL->n + j]);
    }

    for (int k = j + 1; k < SL->n; k++) // After the current target (matrix number).
    {
      alpha += fabs(SL->A[line * SL->n + j]);
    }

    alpha /= SL->A[line * SL->n + j]; // Division by current term A(line, j).

    if (lastAlpha < alpha)
    {
      jAlpha = j;
    }
  }

  return jAlpha;
}

real_t calculaErro(real_t *xAtual, real_t *xAnterior, int n)
{
  real_t erroMaximo = xAtual[0] - xAnterior[0];
  real_t erro;

  for (int i = 1; i < n; i++)
  {
    erro = xAtual[i] - xAnterior[i];

    if (erro > erroMaximo)
    {
      erroMaximo = erro;
    }
  }

  return erroMaximo;
}

/**
 * @brief Método de Gauss-Jacobi.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param x Ponteiro para o vetor solução.
 * @param erro Menor erro aproximado para encerrar as iterações.
 * @return int Código de erro. Um número positivo indica sucesso e o número de iterações realizadas. Um número negativo indica um erro.
 */
int gaussJacobi(SistLinear_t *SL, real_t *x, real_t erro)
{
  memset(x, 0, SL->n * sizeof(real_t)); // Set the initial solution to zero.

  int iter = 0;

  real_t *tmpX = malloc(SL->n * sizeof(real_t));

  for (int i = 0; i < SL->n; i++) // Sorting the columns (dominant diagonal).
  {
    int jAlpha = encontraMaxAlpha(SL, i);

    if (i != jAlpha)
    {
      trocaPosicao(SL, i, jAlpha);
    }
  }

  do // Calculating the new solution.
  {
    memcpy(tmpX, x, SL->n * sizeof(real_t));

    for (int i = 0; i < SL->n; i++)
    {
      x[i] = SL->b[i];

      for (int j = 0; j < i; j++)
      {
        x[i] -= SL->A[i * SL->n + j] * tmpX[j];
      }

      for (int j = i + 1; j < SL->n; j++)
      {
        x[i] -= SL->A[i * SL->n + j] * tmpX[j];
      }

      x[i] /= SL->A[i * SL->n + i];
    }

		double normaL2 = normaL2Residuo(SL, x);
    erro = calculaErro(x, tmpX, SL->n);
    iter++;

  } while (fabs(erro) >= erro && iter < MAXIT); // Partial column pivot (based on alpha value).

  if (iter == MAXIT)
  {
    fprintf(stderr, MAX_IT_MSGE); // Max allowed iterations.

    return MAX_IT_CODE;
  }

  free(tmpX);

  return 0;
}

/**
 * @brief Function to calculate the maximum beta used in the Gauss Jacobi method.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param line Current line.
 * @return int Beta value.
 */
int encontraMaxBeta(SistLinear_t *SL, int line)
{
  real_t *beta = malloc(SL->n * sizeof(real_t));
  int jBeta = 0;

  for (int j = line; j < SL->n; j++)
  {
    beta[j] = 0;

    for (int k = line; k < j; k++)
    {
      beta[j] += fabs(SL->A[line * SL->n + j]) * beta[k];
    }

    for (int k = j + 1; k < SL->n; k++)
    {
      beta[j] += fabs(SL->A[line * SL->n + j]);
    }

    beta[j] /= SL->A[line * SL->n + j];

    if (beta[jBeta] < beta[j])
    {
      jBeta = j;
    }
  }

  return jBeta;
}

/**
 * @brief Método de Gauss-Seidel.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param x Ponteiro para o vetor solução.
 * @param erro Menor erro aproximado para encerrar as iterações.
 * @return int Código de erro. Um número positivo indica sucesso e o número de iterações realizadas. Um número negativo indica um erro.
 */
int gaussSeidel(SistLinear_t *SL, real_t *x, real_t erro)
{
  memset(x, 0, SL->n * sizeof(real_t));
  int iter = 0;

  for (int i = 0; i < SL->n; i++)
  {
    int jBeta = encontraMaxBeta(SL, i);

    if (i != jBeta)
    {
      trocaPosicao(SL, i, jBeta);
    }
  }

  real_t *tmpX = malloc(SL->n * sizeof(real_t));

  do
  {
    memcpy(tmpX, x, SL->n * sizeof(real_t));

    for (int i = 0; i < SL->n; i++)
    {
      x[i] = SL->b[i];

      for (int j = 0; j < i; j++)
      {
        x[i] -= SL->A[i * SL->n + j] * x[j];
      }

      for (int j = i + 1; j < SL->n; j++)
      {
        x[i] -= SL->A[i * SL->n + j] * x[j];
      }

      x[i] /= SL->A[i * SL->n + i];
    }

    double normaL2 = normaL2Residuo(SL, x);
    erro = calculaErro(x, tmpX, SL->n);
    iter++;
  } while (fabs(erro) >= erro && iter < MAXIT);

  if (iter == MAXIT)
  {
    fprintf(stderr, MAX_IT_MSGE); // Max allowed iterations.

    return MAX_IT_CODE;
  }

  free(tmpX);

  return 0;
}

// Alocaçao de memória
SistLinear_t *alocaSistLinear(unsigned int tam)
{
  SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));
  if (SL)
  {
    SL->A = (real_t *)malloc(tam * tam * sizeof(real_t));
    SL->b = (real_t *)malloc(tam * sizeof(real_t));

    if (!(SL->A) || !(SL->b))
      liberaSistLinear(SL);
  }

  SL->n = tam;

  return (SL);
}

// Liberacao de memória
void liberaSistLinear(SistLinear_t *SL)
{
  free(SL->A);
  free(SL->b);
  free(SL);
}

/**
 * @brief Cria coeficientes e termos independentes do SL.
 *
 * @param SL Ponteiro para o sistema linear.
 * @param tipo Tipo de sistema linear a ser criado. Pode ser: comSolucao, eqNula, eqProporcional, eqCombLinear e hilbert.
 * @param coef_max Maior valor para coeficientes e termos independentes.
 */
void inicializaSistLinear(SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max)
{
  unsigned int tam = SL->n;
  // Para gerar valores no intervalo [0, coef_max].
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);

  // Inicializa vetor b.
  for (unsigned int i = 0; i < tam; ++i)
  {
    SL->b[i] = (real_t)rand() * invRandMax;
  }

  if (tipo == hilbert)
  {
    for (unsigned int i = 0; i < tam; ++i)
    {
      for (unsigned int j = 0; j < tam; ++j)
      {
        SL->A[i * tam + j] = 1.0 / (real_t)(i + j + 1);
      }
    }
  }
  else
  { // Inicializa sistema normal e depois altera.
    // Inicializa a matriz A.
    for (unsigned int i = 0; i < tam; ++i)
    {
      for (unsigned int j = 0; j < tam; ++j)
      {
        SL->A[i * tam + j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula)
    {
      // Sorteia eq a ser "nula".
      unsigned int nula = rand() % tam;
      for (unsigned int j = 0; j < tam; ++j)
      {
        SL->A[nula * tam + j] = 0.0;
      }
      SL->b[nula] = 0.0;
    }
    else if (tipo == eqProporcional)
    {
      // Sorteia eq a ser "proporcional" e valor.
      unsigned int propDst = rand() % tam;
      unsigned int propSrc = (propDst + 1) % tam;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j = 0; j < tam; ++j)
      {
        SL->A[propDst * tam + j] = SL->A[propSrc * tam + j] * mult;
      }
      SL->b[propDst] = SL->b[propSrc] * mult;
    }
    else if (tipo == eqCombLinear)
    {
      // Sorteia eq a ser "combLinear".
      unsigned int combDst = rand() % tam;
      unsigned int combSrc1 = (combDst + 1) % tam;
      unsigned int combSrc2 = (combDst + 2) % tam;
      for (unsigned int j = 0; j < tam; ++j)
      {
        SL->A[combDst * tam + j] = SL->A[combSrc1 * tam + j] + SL->A[combSrc2 * tam + j];
      }
      SL->b[combDst] = SL->b[combSrc1] + SL->b[combSrc2];
    }
    else if (tipo == diagDominante)
    {
      // Aumenta o expoente dos termos da diagonal principal.
      for (unsigned int i = 0; i < tam; ++i)
      {
        SL->A[i * tam + i] *= (real_t)tam;
      }
    }
  }
}

SistLinear_t *lerSistLinear()
{
  unsigned int n;
  SistLinear_t *SL;

  scanf("%d", &n);

  SL = alocaSistLinear(n);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      scanf("%g", &SL->A[i * n + j]);

  for (int i = 0; i < n; ++i)
    scanf("%g", &SL->b[i]);

  return SL;
}

void prnSistLinear(SistLinear_t *SL)
{
  int n = SL->n;

  for (int i = 0; i < n; ++i)
  {
    printf("\n");
    for (int j = 0; j < n; ++j)
      printf("%10g", SL->A[i * n + j]);
    printf("   |   %g", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor(real_t *v, unsigned int n)
{
  int i;

  printf("\n");
  for (i = 0; i < n; ++i)
    printf("%10g ", v[i]);
  printf("\n\n");
}
