#include "Vec.h"

/***** FUNCTION MEMBERS *****/

Vec::Vec (int n, double v){
  N = n;
  entries = new double [N];
  for (int i = 0; i < N; i++)
    entries[i] = v;
} // constructor with num. el. and value

Vec::Vec (int n, double * vector){
  N = n;
  entries = new double [n];

  for(int i = 0; i < n; i++)
    entries[i] = vector[i];
} // constructor with num. el. and ptr. array

Vec::Vec (vector<double> vector){
  N = vector.size();
  entries = new double [N];
  for (int i = 0; i < N; i++)
    entries[i] = vector[i];
} // constructor with vector<double>

Vec::~Vec() {
  delete[] entries;
} // destructor

Vec::Vec (const Vec & aceita_e_copia){
  N = aceita_e_copia.N;
  entries = new double [aceita_e_copia.N];

  for(int i = 0; i < N; i++)
    entries[i] = aceita_e_copia.entries[i];
} // copy constructor...

void Vec::SetEntries (int n, double *array){
  if (N == n){
    for (int i = 0; i < n; i++)
      entries[i] = array[i];
  }

  if (N != n){
    N = n;
    delete[] entries;
    entries = new double [n];
    for (int i = 0; i < n; i++)
      entries[i] = array[i];
  }
} // set entries

void Vec::Print(){
  cout << "| ";
  for(int i = 0; i < N - 1; i++){
    cout << setw(10) << entries[i];
  }
  cout << setw(10) << entries[N-1] << " |" << endl;
} // print the vector content

void Vec::swap (int n, int m){
  double c = entries[n];

  entries[n] = entries[m];
  entries[m] = c;
} // swap

int Vec::size () const{
  return N;
} // return size of vector

double Vec::dot (const Vec & obj){
  if (N != obj.N)
    return 0;

  double dot = 0;
  for (int i = 0; i < N; i++)
    dot = dot + obj.entries[i] * entries[i];

  return dot;
} // scalar product with another Vec

Vec & Vec::operator= (const Vec &twins){
  if (this != &twins){
    if (N != twins.N){
      N = twins.N;
      delete[] entries;
      entries = new double[twins.N];

      for( int i = 0; i < N; i++)
        entries[i] = twins.entries[i];
    }

    if (N == twins.N){
      for( int i = 0; i < N; i++)
        entries[i] = twins.entries[i];
    }
  }
  return *this;
} // operator=

Vec Vec::operator+ (const Vec & obj){
  if (N != obj.N){
    cout << __PRETTY_FUNCTION__ << " : Vec's must have the same Dimension (" << N << " != " << obj.N << ")" <<  endl;
    return Vec(N); // return Vec 0 with N elements
  }

  double* sum = new double [N];
  for (int i = 0; i < N; i++)
    sum[i] = entries[i] + obj.entries[i];

  Vec result;
  result.SetEntries(N, sum);
  delete[] sum;

  return result;
} // operator+

Vec & Vec::operator+= (const Vec &soma){

  if (N != soma.N){
    cout << __PRETTY_FUNCTION__ << " : Vec's must have the same Dimension (" << N << " != " << soma.N << ")" <<  endl;
    return *this; // Keep Vec Unmodified
  }

  for(int i = 0; i < N; i++)
    entries[i] += soma.entries[i];

  return *this;
} // operator+=

Vec Vec::operator- (const Vec & obj){
  if (N != obj.N){
    cout << __PRETTY_FUNCTION__ << " : Vec's must have the same Dimension (" << N << " != " << obj.N << ")" <<  endl;
    return Vec(N); // return Default Vec
  }

  double* sub = new double [N];
  for (int i = 0; i < N; i++)
    sub[i] = entries[i] - obj.entries[i];

  Vec result;
  result.SetEntries(N, sub);
  delete[] sub;

  return result;
} // operator-

Vec & Vec::operator-= (const Vec &subtracao){

  if (N != subtracao.N){
    cout << __PRETTY_FUNCTION__ << " : Vec's must have the same Dimension (" << N << " != " << subtracao.N << ")" << endl;
    return *this; // Keep Vec Unmodified
  }

  for(int i = 0; i < N; i++)
    entries[i] += subtracao.entries[i];

  return *this;
} // operator-=

double & Vec::operator[] (int x){
  if (x < 0 || x >= N){
    cout << __PRETTY_FUNCTION__ << " : Dimension Error (" << x << " >= " << N << ")" << endl;
    return entries[N - 1]; // if Err -> return Last Element of Vec
  }
  return entries[x];
} // operator[]

double Vec::operator[] (int x) const{
  if (x < 0 || x >= N){
    cout << __PRETTY_FUNCTION__ << " : O elemento " << x << " não existe." << " O valor admitido foi entries[" << N - 1 << "]." << endl;
    return entries[N - 1];
  }

  return entries[x];
} // operator[] const

Vec Vec::operator* (const Vec & obj) const{
  if (N != obj.N){
    cout << __PRETTY_FUNCTION__ << " : Vec's must have the same Dimension (" << N << " != " << obj.N << ")" <<  endl;
    return Vec(N); // return Vec 0 with N elements
  }

  double* prod = new double [N];
  for (int i = 0; i < N; i++)
    prod[i] = entries[i] * obj.entries[i];

  Vec result;
  result.SetEntries(N, prod);
  delete[] prod;

  return result;
} // operator*

Vec Vec::operator* (const double & scalar) const{ // const function !!
  double* array_auxiliar = new double [N];
  for(int i = 0; i < N ; i++){
    array_auxiliar[i] = entries[i] * scalar;
  }

  Vec result;
  result.SetEntries(N, array_auxiliar);
  delete[] array_auxiliar;

  return result;
} // multiplication by scalar (Vec LEFT)

Vec operator* (const double & scalar, const Vec & V){
  Vec copy_V(V);
  return (copy_V * scalar);
} // multipliation by scalar (Vec RIGHT) -- friend function

Vec Vec::operator+ () const {
  return Vec(*this);
} // operator+ (unary)

Vec Vec::operator- () const {
  vector<double> aux;
  for (int i = 0; i < N; i++){
    aux.push_back(-entries[i]);
  }
  return Vec(aux);
} // operator- (unary)
