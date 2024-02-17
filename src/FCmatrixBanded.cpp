#include "FCmatrixBanded.h"

/************ Constructors ****************/

FCmatrixBanded::FCmatrixBanded(){
    M.push_back(Vec());
    M.push_back(Vec(2));
    M.push_back(Vec());
    classname = "FCmatrixBanded";
} //default constructor

FCmatrixBanded::FCmatrixBanded(double** fM, int fm){
  vector<double> aux;
  int i;

  for (i = 0; i < fm - 1 ; i++){
    aux.push_back(fM[i][i+1]);
  }
  M.push_back(Vec(aux));
  aux.clear();

  for (i = 0; i < fm; i++){
    aux.push_back(fM[i][i]);
  }
  M.push_back(Vec(aux));
  aux.clear();

  for (i = 0; i < fm - 1; i++){
    aux.push_back(fM[i+1][i]);
  }
  M.push_back(Vec(aux));

  classname = "FCmatrixBanded";
}

FCmatrixBanded::FCmatrixBanded(double* fM, int fm){
  vector<double> aux;
  int i;

  for ( i = 0; i < fm - 1; i++){
    aux.push_back(fM[i]);
  }

  M.push_back(Vec(aux));
  aux.clear();

  while (i < 2*fm - 1){
    aux.push_back(fM[i]);
    i++;
  }

  M.push_back(Vec(aux));
  aux.clear();

  while (i < 3*fm - 2){
    aux.push_back(fM[i]);
    i++;
  }

  M.push_back(Vec(aux));

  classname = "FCmatrixBanded";
} //matrix fm x fn

FCmatrixBanded::FCmatrixBanded(vector<Vec> V): FCmatrix(V){
  classname = "FCmatrixBanded";
} // matrix fm x fn given as vector of Vec

FCmatrixBanded::FCmatrixBanded (const FCmatrixBanded & copy ){
    M = copy.GetM();
    classname = "FCmatrixBanded";
} // copy constructor

void FCmatrixBanded::Print() const{
  FCmatrixFull conversao;
  conversao = convert(*this);
  conversao.Print();
} // Print de M

/****************** Operators *****************/

FCmatrixBanded FCmatrixBanded::operator=(const FCmatrixBanded & copy){
  if (this != &copy){
    classname = copy.GetClassname();
    M = copy.GetM();
  }
  return (*this);
} // it doesn't work with FCmatrix as argument ...


FCmatrixBanded FCmatrixBanded::operator+(const FCmatrixBanded & soma) const{
    vector<Vec> aux;
    if(Get_nRows() != soma.Get_nRows()){
      cout << __PRETTY_FUNCTION__ << " : Err." << endl;
      return FCmatrixBanded();
    }

    for (int i=0; i<3; i++){
      aux.push_back(GetDiagonal(i) + soma.GetDiagonal(i));
    }
    return FCmatrixBanded(aux);
}


FCmatrixBanded FCmatrixBanded::operator-(const FCmatrixBanded & sub) const{ // sub 2 matrices of any kind
  vector<Vec> aux;

  if(Get_nRows() != sub.Get_nRows()){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return FCmatrixBanded();
  }

  for (int i=0; i<3; i++){
    aux.push_back(GetDiagonal(i) - sub.GetDiagonal(i));
  }

  return FCmatrixBanded(aux);
}


FCmatrixBanded FCmatrixBanded::operator*(double lambda) const{
  vector<Vec> aux;

  for (int i=0; i<3; i++){
    aux.push_back(GetDiagonal(i) * lambda);
  }

  return FCmatrixBanded(aux);
} // mul matrix of any kind by scalar (Matrix LEFT)


FCmatrixBanded operator* (double lambda, const FCmatrixBanded & FCB){
  FCmatrixBanded result = FCB * lambda;
  return result;
} // mul matrix of any kind by scalar (Matrix RIGHT) --- NON MEMBER FUNCTION


Vec FCmatrixBanded::operator* (const Vec & col_vec) const{
  FCmatrixFull Full_M = convert((*this));
  return (Full_M * col_vec);
} // mul matrix by Vec


int FCmatrixBanded:: Get_nRows() const{
  if (GetM().size() == 1)
    return 1;

  return M[1].size();
}


int FCmatrixBanded:: Get_nCols() const{
  if (GetM().size() == 1)
    return 1;

  return M[1].size();
}


Vec FCmatrixBanded:: GetRow(int i) const{
  return convert(*this).GetRow(i);
}


Vec FCmatrixBanded:: GetCol(int i) const{
  return convert(*this).GetCol(i);
}


Vec FCmatrixBanded:: GetDiagonal(int i) const{
  /* i = 0 : GetUpperDiagonal()
     i = 1 : GetMainDiagonal()
     i = 2 : GetLowerDiagonal() */
  if (Get_nRows() == 1)
    return M[0];

  return M[i];
}


Vec& FCmatrixBanded:: operator[] (int i){
  if (i > 2 || i < 0){
    cout << __PRETTY_FUNCTION__ << " : Err." << endl;
    return M[2];
  }
  return M[i];
}


FCmatrixFull convert (const FCmatrixBanded origin){
  int n = origin.Get_nRows();
  vector<Vec> aux;
  for (int i = 0; i < n ; i++){
    aux.push_back(Vec(n));
  }
  // Diagonal
  for (int i = 0; i < n; i++){
    aux[i][i]= origin.GetDiagonal(1)[i];
  }

  // Diagonal Superior
  for (int i=0; i < n - 1; i++){
    aux[i][i+1] = origin.GetDiagonal(0)[i];
  }

  // Diagonal Inferior
  for (int i=0; i < n - 1; i++){
    aux[i+1][i] = origin.GetDiagonal(2)[i];
  }
  return FCmatrixFull(aux);
} // NON MEMBER FUNCTION (friend conversion)
