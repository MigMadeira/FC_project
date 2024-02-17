#include "FCmatrix.h"

FCmatrix::FCmatrix(){
  classname = "FCmatrix";
} // default constructor

FCmatrix::FCmatrix(double** fM, int fm, int fn){
  for (int i = 0; i < fm; i++){
    M.push_back(Vec(fn, fM[i]));
  }
  classname = "FCmatrix";
} // matrix fm xfn given from pointer of pointers

FCmatrix::FCmatrix (double* fM, int fm, int fn){
  vector<double> vector;

  for (int i = 0; i <= fn*fm; i++){
    if (i % fn == 0 && i != 0){
      M.push_back(Vec(vector));
      vector.clear();
    }
    if (i != fn*fm)
      vector.push_back(fM[i]);
  }
} // matrix fm x fn given as single pointer (what length ?!) -- length : fm*gn

FCmatrix::FCmatrix(vector<Vec> V){
  M = V;
  classname = "FCmatrix";
} // matrix fm x fn given as vector of Vec

void FCmatrix::Print() const{
  for (int i = 0; i < M.size(); i++){
    (GetRow(i)).Print();
  }
} // Print de M

vector<Vec> FCmatrix::GetM() const{
  return M;
} // access to matrix M

string FCmatrix::GetClassname() const{
  return classname;
} // access to classname
