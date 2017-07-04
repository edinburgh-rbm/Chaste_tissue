#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "TextMatrix.hpp"

using namespace std;

double **read_matrix(const char* filename, int& rows, int&cols) {
  ifstream fin;
 
  fin.open (filename);
  if (!fin) return NULL;
 
  fin >> rows >> cols;
 
  double **a = new double *[rows]; 
 
  for (int i = 0; i < rows; i++) {
    a[i] = new double[cols];
  }
  
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      fin >> a[i][j];
    }
  }

  return a;
}
