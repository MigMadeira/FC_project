#ifndef H_EQSOLVER_H
#define H_EQSOLVER_H
#include "FCmatrixFull.h"
#include "FCmatrixBanded.h"

class EqSolver {
  public:
    EqSolver();
    EqSolver(const FCmatrixFull &, const Vec &); // matriz M e vector de constantes
    EqSolver(const FCmatrixBanded &, const Vec &); // matriz tridiagonal M e vector de constantes

    // set
    void SetConstants (const Vec &);
    void SetMatrix (const FCmatrixFull &);
    void SetMatrix (const FCmatrixBanded &);

    // Solve Functions
    Vec GaussEliminationSolver();
    Vec LUdecompositionSolver_Doolittle();
    Vec TridiagonalSolver();
    // Vec JacobiSolver(double tol = 1.E-6);

  private:
    // return triangular matrix and changed vector of constants
    void GaussElimination(FCmatrixFull &, Vec &);
    // decomposição LU com |L|=1
    FCmatrixFull* LUdecomposition_Doolittle(FCmatrixFull &);
    // void LUdecomposition(FCmatrixFull &, vector<int> & index);
    Vec TridiagonalThomas(FCmatrixBanded &, Vec &);

    FCmatrixFull M; //matriz de coeffs
    FCmatrixBanded Band; //objecto com bandas
    Vec b; //vector de constantes

    bool isSolved;
    bool isLUSolved;
    bool isMatrix; //full matrix set ?
    bool isBand; //banded matrix set ?
    bool isB; //vector b set ?
};

#endif
