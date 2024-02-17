#include "FCmatrixFull.h"

FCmatrixFull::FCmatrixFull(){
  M.push_back(Vec()); // default Vec()
  classname = "FCmatrixFull";
} // default constructor

FCmatrixFull::FCmatrixFull(double** fM, int fm, int fn): FCmatrix(fM, fm, fn){
  classname = "FCmatrixFull";
} // matrix fm xfngiven from pointer of pointers

FCmatrixFull::FCmatrixFull(double* fM, int fm, int fn): FCmatrix(fM, fm, fn){
  classname = "FCmatrixFull";
} // matrix fm x fn given as single pointer (what length ?!) -- length : fm*gn

FCmatrixFull::FCmatrixFull(vector<Vec> V): FCmatrix(V){
  classname = "FCmatrixFull";
} // matrix fm x fn given as vector of Vec

FCmatrixFull::FCmatrixFull (const FCmatrixFull & origin){
  classname = origin.GetClassname();
  M = origin.GetM();
} // copy constructor

/******************* OPERATORS ********************/

FCmatrixFull FCmatrixFull::operator=(const FCmatrix & copy){
  if (this != &copy){
    classname = copy.GetClassname();
    M = copy.GetM();
  }
  return (*this);
} // equal 2 matrices of any kind (copy assignement)

FCmatrixFull FCmatrixFull::operator+(const FCmatrix & matrix2) const{
  vector<Vec> aux;
  if (Get_nRows() != matrix2.Get_nRows() || Get_nCols() != matrix2.Get_nCols()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return FCmatrixFull();
  }

  for(int i = 0; i < Get_nRows(); i++){
    aux.push_back(GetRow(i) + matrix2.GetRow(i));
  }
  return FCmatrixFull (aux);
} // add 2 matrices of any kind

FCmatrixFull FCmatrixFull::operator-(const FCmatrix & obj) const{
  if (obj.Get_nRows() != Get_nRows() || obj.Get_nCols() != Get_nCols()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return FCmatrixFull();
  }

  vector<Vec> M_sub;
  for (int i = 0; i < obj.Get_nRows(); i++){
    M_sub.push_back(GetRow(i) - obj.GetRow(i));
  }
  return FCmatrixFull(M_sub);
} // sub 2 matrices of any kind

FCmatrixFull FCmatrixFull::operator*(const FCmatrix & obj) const{
  if (Get_nRows() != obj.Get_nCols()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return FCmatrixFull();
  }

  vector<Vec> M_prod;
  vector<double> aux;

  for (int i = 0; i < Get_nRows(); i++){
    for (int j = 0; j < obj.Get_nCols(); j++){
      aux.push_back((GetRow(i)).dot(obj.GetCol(j)));
    }

    M_prod.push_back(Vec(aux));
    aux.clear();
  }
  return FCmatrixFull(M_prod);
} // mul 2 matrices of any kind

Vec FCmatrixFull::operator*(const Vec & b) const{
  if (Get_nCols() != b.size()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return b; // return 'b' if dimensions don't match
  }
  vector<double> result;
  for (int i = 0; i < M.size(); i++){
    result.push_back(GetRow(i).dot(b));
  }
  return Vec(result);
} // mul matrix by Vec

FCmatrixFull FCmatrixFull::operator*(double lambda) const{
  vector<Vec> result;

  for (int i = 0; i < M.size(); i++){
    result.push_back(lambda * GetRow(i));
  }

  return FCmatrixFull(result);
} // mul matrix of any kind by scalar (Matrix LEFT)

FCmatrixFull operator* (double lambda, const FCmatrixFull & FMF){
  FCmatrixFull result(FMF * lambda);
  return result;
} // mul matrix of any kind by scalar (Matrix RIGHT) --- NON MEMBER FUNCTION (friend)

/**************** VIRTUAL INHERITED ****************/

int FCmatrixFull::Get_nRows() const{
  return M.size();
} // return number of rows

int FCmatrixFull::Get_nCols() const {
  return M[0].size();
} // return number of columns

Vec FCmatrixFull::GetRow(int i) const{
  return M[i];
} // return Row(index i)

Vec FCmatrixFull::GetCol(int i) const{
  vector<double> aux;
  for (int j = 0; j < M.size(); j++){
    aux.push_back(M[j][i]);
  }
  return Vec(aux);
} // return Col(index i)

double FCmatrixFull::Determinant() const{

} // determinant (soon...)


Vec& FCmatrixFull::operator[] (int i){
  if (i < 0 || i >= M.size())
    return M[M.size() - 1]; // by default it returns the last Vec of the vector<Vec>

  return M[i];
}; // operator[]

int FCmatrixFull::GetRowMax (int i) const{
  if (i >= M.size()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return 0; // default return 0
  }

  int max = 0;
  double value = abs(M[i][0]);

  for (int j = 0; j < M[i].size(); j++){
    if (value < abs(M[i][j])){
      value = abs(M[i][j]);
      max = j;
    }
  }
  return max;
} // in row-i, return column-index of max element (in absolute value)

int FCmatrixFull::GetColMax(int j) const{
  // Se a matriz não for quadrada e o número de linhas for menor
  // que o número de colunas, apenas é possível determinar o index-ColMax para, no máximo, j < M[0].size()

  if (j >= M.size() || j >= M[0].size()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return 0; // default return 0
  }

  int k = j; // Encontrar a primeira linha não nula
  while (M[k][GetRowMax(k)] == 0 && k < Get_nRows()){
    k++;
  }

  if (k == (M.size() - 1)) // Se as linhas forem todas nulas retornar a última linha
    return k;

  int max = k;
  double scaled_value = abs(M[k][j] / M[k][GetRowMax(k)]);

  for (int i = k + 1; i < M.size(); i++){
    if (GetRow(i)[GetRowMax(i)] == 0){
      continue;
    }

    if (scaled_value < abs(M[i][j] / M[i][GetRowMax(i)])){
      max = i;
      scaled_value = abs(M[i][j] / M[i][GetRowMax(i)]);
      continue;
    }

    /*
    if (scaled_value == abs(M[i][j] / M[i][GetRowMax(i)])){ // Non-scaled Partial Pivoting (if necessary ...)
                                                           //  Check example array_B on Test_EQ.cpp
      if (abs(M[i][j]) > abs(M[max][j])){
        max = i;
        // scaled_value = abs(M[i][k] / M[i][GetRowMax(i)]);
        continue;
      }
    }
    */
  }

  /* For non-scaled/relative Partial Pivoting :
  int max = j;
  double scaled_value = abs(M[j][j]);
  for (int i = j; i < M.size(); i++){
    if (scaled_value < abs(M[i][j])){
      max = i;
      scaled_value = abs(M[i][j]);
    }
  } */

  return max;
} // in column-j, return row-index (>=j) for which relative amplitude of Mij on the row is highest.

void FCmatrixFull::swapRows(int i, int j){
  Vec aux;
  aux = M[i];
  M[i] = M[j];
  M[j] = aux;
} // swap Rows index i and j
