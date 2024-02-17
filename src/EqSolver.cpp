#include "EqSolver.h"

EqSolver::EqSolver(){
  isMatrix = false;
  isBand = false;
  isB = false;
  isSolved = false;
  isLUSolved = false;
} // default constructor

EqSolver::EqSolver (const FCmatrixFull & FCF, const Vec & Vec){
  M = FCF;
  b = Vec;
  isMatrix = true;
  isBand = false;
  isB = true;
  isSolved = false;
  isLUSolved = false;
} // matriz M e vector de constantes

EqSolver::EqSolver (const FCmatrixBanded & FCB, const Vec & Vec){
  Band = FCB;
  b = Vec;
  isMatrix = false;
  isBand = true;
  isB = true;
  isSolved = false;
  isLUSolved = false;
} // matriztridiagonal M e vector de constantes */

void EqSolver::SetConstants (const Vec & Vec){
  b = Vec;
  if (!(isB)) isB = true;
  isSolved = false;
  isLUSolved = false;
}

void EqSolver::SetMatrix (const FCmatrixFull & FCF){
  M = FCF;
  if (!(isMatrix)) isMatrix = true;
  isSolved = false;
  isLUSolved = false;
}

void EqSolver::SetMatrix (const FCmatrixBanded & FCB){
  Band = FCB;
  if (!(isBand)) isBand = true;
  isSolved = false;
  isLUSolved = false;
}

/******************* DOOLITTLE *********************/

Vec EqSolver::LUdecompositionSolver_Doolittle(){
  if (!(isMatrix) || isSolved || isLUSolved || !(isB)){
    cout << __PRETTY_FUNCTION__ << " : Err. (Check Equation Status)" << endl;
    return Vec();
  }

  if (M.Get_nCols() != M.Get_nRows()){
    cout << __PRETTY_FUNCTION__ << " : Square Matrix Required." << endl;
    return Vec();
  }

  FCmatrixFull* LU = LUdecomposition_Doolittle(M);
  /* LU[0].Print();
  LU[1].Print(); */

  Vec y(b.size());
  // Ly = b
  for (int i = 0; i < LU[0].Get_nRows(); i++)
    y[i] = (b[i] - LU[0][i].dot(y)); // L[i][i] = 1 by definition

  Vec x(b.size());
  // back substitution Ux = y
  for (int i = LU[1].Get_nRows() - 1; i >= 0; i--){
    // check LU[1][i][i] == 0
    if (LU[1][i][i] == 0){
      cout << __PRETTY_FUNCTION__ << " : Can't solve by Doolittle (U[" << i << "][" << i << "] = 0)." << endl;
      delete[] LU;
      isLUSolved = true;
      return x;
    }
    x[i] = (y[i] - LU[1][i].dot(x)) / LU[1][i][i];
  }

  delete[] LU;

  isLUSolved = true;
  return x;
}


FCmatrixFull* EqSolver::LUdecomposition_Doolittle (FCmatrixFull & A){
  int n = A.Get_nRows(); // == Get_nCols() ...

  vector<Vec> vector_L;
  vector<Vec> vector_U;

  for (int i = 0; i < n; i++){
    vector_L.push_back(Vec(n));
    vector_L[i][i] = 1;
    vector_U.push_back(Vec(n));
  }

  FCmatrixFull L(vector_L);
  FCmatrixFull U(vector_U);

  for (int i = 0; i < n; i++){
    for (int j = i; j < n; j++){
      U[i][j] = A[i][j] - L.GetRow(i).dot(U.GetCol(j));
    }

    for (int j = i + 1; j < n; j++){
      if (U[i][i] == 0){
        cout << __PRETTY_FUNCTION__ << " : Skip (U[" << i << "][" << i << "] = 0)." << endl;
        continue;
      }
      L[j][i] = (A[j][i] - L.GetRow(j).dot(U.GetCol(i))) / U[i][i];
    }
  }

  FCmatrixFull* LU = new FCmatrixFull [2];
  LU[0] = L;
  LU[1] = U;

  /* LU[0].Print();
  cout << endl;
  LU[0].Print();
  cout << endl; */

  return LU;
}

/************* GAUSS ELIMINATION ****************/

Vec EqSolver::GaussEliminationSolver (){
  if (!(isMatrix) || isSolved || !(isB)){
    cout << __PRETTY_FUNCTION__ << " : Err. (Check Equation Status)" << endl;
    return Vec();
  }

  GaussElimination(M, b);

  /* M.Print();
  cout << endl;
  b.Print();
  cout << endl; */

  if (isSolved == true || M.Get_nRows() != M.Get_nCols()){
    cout << __PRETTY_FUNCTION__ << " : O Sistema é Impossível ou Indeterminado (Can't solve it ...)" << endl;
    cout << "  -- Gauss Elimination Result --  " << endl << endl;
    M.Print();
    cout << endl;
    b.Print();
    cout << endl;
    isSolved = true; // necessary for non-squared matrices
    return Vec();
  }

  Vec x(b.size());

  for (int i = M.Get_nRows() - 1; i >= 0; i--){
    // check if system is undetermined or impossible
    if (M[i][i] == 0){
      if (M[i][M.GetRowMax(i)] != 0){
        cout << __PRETTY_FUNCTION__ << " : Gauss Elimination (Numeric Err.) : " << endl << endl;
        M.Print();
        cout << endl;
        isSolved = true;
        return Vec();

      }
      if (b[i] == 0){
        cout << __PRETTY_FUNCTION__ << " : O Sistema é Indeterminado." << endl;
        isSolved = true;
        return Vec();
      }
      cout << __PRETTY_FUNCTION__ << " : O Sistema é Impossível." << endl;
      isSolved = true;
      return Vec();
    }
    x[i] = (b[i] - M[i].dot(x)) / M[i][i];
  }

  isSolved = true;

  cout << "  -- Gauss Elimination Result --  " << endl << endl;
  M.Print();
  cout << endl;
  b.Print();
  cout << endl;

  return x;
}


// Private Solver Step
void EqSolver::GaussElimination (FCmatrixFull & A, Vec & b){
  int nrows = A.Get_nRows();
  int ncols = A.Get_nCols();
  double lambda;
  int drow;

  for (int row = 0; row < nrows - 1 && row < ncols; row++){
    drow = A.GetColMax(row);
    /* A.Print();
    cout << endl;
    b.Print();
    cout << endl;

    cout << "row : " << row << endl;
    cout << "drow : " << drow << endl; */

    if (drow == (A.Get_nRows() - 1) && A[drow][A.GetRowMax(drow)] == 0){
      /* se o GetColMax der return de k (ver def.), significa que o sistema tem todas as linhas i <= row
      iguais a 0 */
      break;
    }

    if (A[drow][row] == 0){
      // se, a partir da linha 'row', a coluna 'row' já estiver toda a 'zeros', então o sistema é ind. ou imp.
      isSolved = true;
      break;
    }

    A.swapRows(row, drow);
    b.swap(row, drow);

    for (int i = row + 1; i < nrows; i++){
      lambda = A[i][row] / A[row][row]; // scaling factor
      A[i] = A[i] - A[row] * lambda; // calculate the new row-i
      b[i] = b[i] - b[row] * lambda; // calculate new vector b index-i
    }
  }
}

/************************* TRIDIAGONAL THOMAS *********************/

Vec EqSolver::TridiagonalSolver(){
  if (!(isBand) || isSolved || !(isB)){
    cout << __PRETTY_FUNCTION__ << " : Err. (Check Equation Status)" << endl;
    return Vec();
  }

  return TridiagonalThomas(Band, b);
}

Vec EqSolver::TridiagonalThomas(FCmatrixBanded & M, Vec & b){
  int i = 1;
  Vec x(b.size());
  double lambda;

  // 1x1 "Tridiagonal Matrix" Case !!
  if (M.GetM().size() == 1){
    if (M[0][0] == 0){
      cout << __PRETTY_FUNCTION__ << " : Can't solve by Thomas (Err.)" << endl;
      return x;
    }

    x[0] = b[0] / M[0][0];
    return x;
  }

  // GENERAL CASE
  if (M[1][0] == 0){
    cout << __PRETTY_FUNCTION__ << " : Can't solve by Thomas (Err.)" << endl;
    return x;
  }

  M[0][0] = M[0][0] / M[1][0];
  b[0] = b[0] / M[1][0];

  for (i = 1; i < M.Get_nRows() - 1; i++){
    lambda = M[1][i] - M[0][i-1]*M[2][i-1];
    if (lambda == 0){
      cout << __PRETTY_FUNCTION__ << " : Can't solve by Thomas (Err.)" << endl;
      return x;
    }
    M[0][i] = M[0][i] / lambda;
    b[i] = (b[i] - b[i-1]*M[2][i-1]) / lambda;
  }

  lambda = M[1][i] - M[0][i-1]*M[2][i-1];
  if (lambda == 0){
    cout << __PRETTY_FUNCTION__ << " : Can't solve by Thomas (Err.)" << endl;
    return x;
  }

  b[i] = (b[i] - b[i-1]*M[2][i-1]) / lambda;

  /* M.Print();
  cout << endl;
  b.Print();
  cout << endl; */

  x[x.size()-1] = b[b.size()-1];
  for (int i = M.Get_nRows() - 2; i >= 0; i--){
    x[i] = b[i] - x[i+1]*M[0][i];
  }

  isSolved = true;
  return x;
}
