#include "Functions/functions.h"

#include <armadillo>

#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

int GCT::calculate_eigenvalues_3x3_sym(std::vector<double> &inVector, int startMatrix, std::vector<double> &outVector){

    //Assumes a 3x3 symmetric matrix, taking the top half of the matrix, including the diagonal
    //It also assumes it gets the full vector with all values, not only the 3x3 matrix you want
    //Thus startMatrix is also passed to find the start of the matrix.
    //This function is as dumb as they come, and assumes the theory is correct.

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




std::vector<double> GCT::calculate_vector_rotation(std::vector<double> b_before, std::vector<double> b_after, double &theta){
    // Takes two vectors, returns the axis the rotation is about, as well as the angle of rotation theta

    GCT::normalize_vector(b_before); GCT::normalize_vector(b_after);        // Work with normalized vectors

    theta = std::acos(GCT::vector_dot_product(b_before, b_after));          // Angle of rotation in the plane

    return GCT::vector_cross_product(b_before, b_after);                    // Normal vector to the plane of rotation
}

int GCT::normalize_vector(std::vector<double> &A){
    double amp = GCT::vector_amplitude(A);
    
    for(int i = 0; i < A.size(); i++){
        A[i] = A[i]/amp;
    }

    return 0;
}


int GCT::rotate_vector_in_plane(std::vector<double> &vect, const std::vector<double> &axis, double theta){
    
    double cost = std::cos(theta);
    double cost1 = 1-cost;
    double sint = std::sin(theta);

    arma::vec v(vect);

    arma::mat R = {
        {   cost + axis[0]*axis[0]*cost1,   axis[0]*axis[1]*cost1 - axis[2]*sint,   axis[0]*axis[2]*cost1 + axis[1]*sint    },
        {   axis[1]*axis[0]*cost1 + axis[2]*sint,   cost + axis[1]*axis[1]*cost1,   axis[1]*axis[2]*cost1 - axis[0]*sint    },
        {   axis[2]*axis[1]*cost1 - axis[1]*sint,   axis[2]*axis[1]*cost1 + axis[0]*sint,   cost + axis[2]*axis[2]*cost1    }
    };

    v = R*v;

    vect = { v[0], v[1], v[2] };

    return 0;
}

int GCT::scalar_dot_vector(double s, std::vector<double> &v){
  for(int i = 0; i < v.size(); i++){
    v[i] = s*v[i];
  }
    return 0;
}

double GCT::frexp10(double arg, int * exp)
{
   *exp = (arg == 0) ? 0 : 1 + (int)std::floor(std::log10(std::fabs(arg) ) );
   return arg * std::pow(10 , -(*exp));    
}




std::string GCT::generate_unique_filename_positions(Bfield &bfield, Particle &particle, int procID){

    std::string filename;

    int expntn;

    frexp10(particle.E, &expntn);

    filename = "Data/Ee" + std::to_string(expntn-1)  + "_LM" + std::to_string(static_cast<int>(bfield.lambda_max))
                         + "/r_rank" + std::to_string(procID) + ".dat";

    return filename;
}
std::string GCT::generate_unique_filename_eigenvalues(double E, double L_max, int procID){

    std::string filename;

    int expntn;

    GCT::frexp10(E, &expntn);

    filename = "Data/Ee" + std::to_string(expntn-1)  + "_LM" + std::to_string(static_cast<int>(L_max))
                         + "/eigenvalues.dat";

    return filename;
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



void GCT::print_Dij(std::vector<double> &inVector, int startMatrix){

  std::cout << std::setprecision(5);
  std::cout << std::endl;

  std::cout << inVector[startMatrix + 0] << ' ' << inVector[startMatrix + 1] << ' ' << inVector[startMatrix + 2] << std::endl;
  std::cout << inVector[startMatrix + 1] << ' ' << inVector[startMatrix + 3] << ' ' << inVector[startMatrix + 4] << std::endl;
  std::cout << inVector[startMatrix + 2] << ' ' << inVector[startMatrix + 4] << ' ' << inVector[startMatrix + 5] << std::endl;
  
  std::cout << std::endl;

}


std::ostream& operator<<(std::ostream& os, const std::vector<double> &v){
  os << v[0] << ' ' << v[1] << ' ' << v[2];
  return os;
}