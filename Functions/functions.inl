#include "Functions/functions.h"

#include <armadillo>

#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;


template<class T, class L>
int GCT::scalar_dot_vector(T s, std::array<L, 3> &v){
  for(int i = 0; i < v.size(); i++){
    v[i] = s*v[i];
  }
    return 0;
}


template<class T, class L>
std::array<T, 3> GCT::scalar_dot_vector_r(L s, const std::array<T, 3> &v){
  std::array<T, 3> b;
  for(int i = 0; i < v.size(); i++){
    b[i] = s*v[i];
  }
    return b;
}

template<class T>
int GCT::normalize_vector(std::array<T, 3> &A){
    double amp = GCT::vector_amplitude(A);
    
    if(amp == 0) return 1;

    for(int i = 0; i < A.size(); i++){
        A[i] = A[i]/amp;
    }

    return 0;
}