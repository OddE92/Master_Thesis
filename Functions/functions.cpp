#include "Functions/functions.h"

#include <armadillo>

#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

int GCT::calculate_eigenvalues_3x3_sym(std::vector<double> &inVector, int startMatrix, std::vector<double> &outVector){

    // Assumes a 3x3 symmetric matrix, taking the top half of the matrix, including the diagonal
    // It also assumes it gets the full vector with all values, not only the 3x3 matrix you want
    // Thus startMatrix is also passed to find the start of the matrix.
    // This function is as dumb as they come, and assumes the theory is correct.

    // This is used together with the eigenvectors, and is supposed to take vectors (not array<double, 3>)

    arma::mat A(3, 3, arma::fill::zeros);

    A(0,0) = inVector[startMatrix + 0];   A(1,0) = inVector[startMatrix + 1];   A(2,0) = inVector[startMatrix + 2];
    A(0,1) = inVector[startMatrix + 1];   A(1,1) = inVector[startMatrix + 3];   A(2,1) = inVector[startMatrix + 4];
    A(0,2) = inVector[startMatrix + 2];   A(1,2) = inVector[startMatrix + 4];   A(2,2) = inVector[startMatrix + 5];

    arma::vec b(3);

    arma::eig_sym(b, A);

    for(int i = 0; i < 3; i++){
        outVector[i] = b(i);
    }
    
    return 0;
}

template<class T>
double GCT::frexp10(T arg, int * exp)                  // Gets the exponent of arg and stores it in exp
{
   double argg = static_cast<double>(arg);
   *exp = (argg == 0) ? 0 : 1 + (int)std::floor(std::log10(std::fabs(argg) ) );
   return argg * std::pow(10 , -(*exp));    
}



    // Generates a unique filename for storing the positions. Always ends in ~/r_rank[procID].dat
std::string GCT::generate_unique_filename_positions(Bfield &bfield, Particle &particle, int procID){

    std::string filename;

    int expntn;

    frexp10(particle.E, &expntn);

    filename = "Data/Ee" + std::to_string(expntn-1)  + "_LM" + std::to_string(static_cast<int>(bfield.lambda_max))
                         + "/r_rank" + std::to_string(procID) + ".dat";

    return filename;
}

    // Generates a unique filename for storing the eigenvalues. Always ends in ~/eigenvalues.dat
std::string GCT::generate_unique_filename_eigenvalues(double E, double L_max, int procID){

    std::string filename;

    int expntn;

    GCT::frexp10(E, &expntn);

    filename = "Data/Ee" + std::to_string(expntn-1)  + "_LM" + std::to_string(static_cast<int>(L_max))
                         + "/eigenvalues.dat";

    return filename;
}

std::string GCT::generate_unique_filename_compare(double E, int Bpercent, int procID){

    int expntn;

    GCT::frexp10(E, &expntn);

    std::string filename = "Data/Compare/Ee" + std::to_string(expntn-1) + "_Bpercent" + std::to_string(Bpercent);
    
    if(procID == 0) filename += "_GCT";

    return filename + ".dat";
}



int GCT::create_directory_to_file(std::string filename){
    // Expects a path on the form "firstDirectory/secondDirectory/.../file"
    // Checks and makes all directories up until /file. The file is not created.

    std::vector<std::string> directories;
    std::string temp;

    std::stringstream tokenize(filename);

    while(getline(tokenize, temp, '/')){
        directories.push_back(temp);
    }

    for(int i = 0; i < directories.size(); i++){

        if(i==0){
            temp = directories[i];
        } 
        else if(i < directories.size()-1){
            temp = temp + '/' + directories[i];
        }

    }

    std::filesystem::path p = temp;

    std::filesystem::create_directories(p);

    return 0;
}


    // Print function for testing
void GCT::print_Dij(std::vector<double> &inVector, int startMatrix){

  std::cout << std::setprecision(5);
  std::cout << std::endl;

  std::cout << inVector[startMatrix + 0] << ' ' << inVector[startMatrix + 1] << ' ' << inVector[startMatrix + 2] << std::endl;
  std::cout << inVector[startMatrix + 1] << ' ' << inVector[startMatrix + 3] << ' ' << inVector[startMatrix + 4] << std::endl;
  std::cout << inVector[startMatrix + 2] << ' ' << inVector[startMatrix + 4] << ' ' << inVector[startMatrix + 5] << std::endl;
  
  std::cout << std::endl;

}

double GCT::calculate_R_larmor(const mps &v_perp, const eV &E, const microGauss &B_amp){
    return (1.0810076e-15 * (v_perp / GCT::c) * (E / B_amp) );
}

    // Overload << for simple vector-printing while debugging.
std::ostream& operator<<(std::ostream& os, const std::vector<double> &v){
  os << v[0] << ' ' << v[1] << ' ' << v[2];
  return os;
}
std::ostream& operator<<(std::ostream& os, const std::array<double, 3> &v){
  os << v[0] << ' ' << v[1] << ' ' << v[2];
  return os;
}