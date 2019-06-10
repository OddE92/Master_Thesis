#ifndef FREE_FUNCTIONS
#define FREE_FUNCTIONS

#include "Bfield/class_bfield.h"
#include "Particle/class_particle.h"
#include "Initializer/initializer.h"

#include <boost/mpi.hpp>

#include <vector>
#include <string>

namespace GCT{        // Put functions in a namespace for safety-reason


int calculate_eigenvalues_3x3_sym(std::vector<double> &inVector, int startMatrix, std::vector<double> &outVector);


double frexp10(double arg_in, int * exponent_out);        // takes a double and stores the exponent in exponent_out

double calculate_R_larmor(const double v_perp, const double E, const double B_amp);

/********** INLINED VECTOR FUNCTIONS **********/


inline double vector_amplitude(const std::array<double, 3> &vect){
  return(std::sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]));
}

inline std::array<double, 3> vector_cross_product(const std::array<double, 3> &A, const std::array<double, 3> &B){
  return { A[1]*B[2] - A[2]*B[1], A[2]*B[0] - A[0]*B[2], A[0]*B[1] - A[1]*B[0] };
}

inline double vector_dot_product(const std::array<double, 3> &A, const std::array<double, 3> &B){
  return (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
}


/***********************************************/

int scalar_dot_vector(double s, std::array<double, 3> &v);
std::array<double, 3> scalar_dot_vector_r(double s, const std::array<double, 3> &v);



    // Filename generators
std::string generate_unique_filename_positions(Bfield &bfield, Particle &particle, int procID);
std::string generate_unique_filename_eigenvalues(double E, double L_max, int procID);
std::string generate_unique_filename_compare(double E, int B_0, int procID);


    // Directory generator
int create_directory_to_file(std::string filename);


    // Vector functions

std::array<double, 3> calculate_vector_rotation(std::array<double, 3> b_before, std::array<double, 3> b_after, double &theta);
int rotate_vector_in_plane(std::array<double, 3> &vector, const std::array<double, 3> &axis, double theta);
int normalize_vector(std::array<double, 3> &A);


    // For printing diffusion tensor while debugging
void print_Dij(std::vector<double> &inVector, int startMatrix);

}// End namespace GCT


/******************* GENERATE_SAMPLES-FUNCTIONS ********************/

    struct Eigenvectors{
      std::vector< std::vector<double> > all_Da_ij;               // Holds all D^(a)_ij sent from each instance of the MF
      std::vector<double> Da_ij;                                  // Holds one D^(a)_ij, used to store recieved data into all_Da_ij
      std::vector<double> D_average;                              // Holds the average values of D_ij
     
        //   values_by_instace_and_time holds the eigenvalues for each MF-instance, at all the recorded points in time
        //   values_at_time holds the eigenvalue summed over all instances, at all the recorded points in time
        //   current_value is a placerholder used when calculating the eigenvalues

        std::vector< std::vector< std::vector<double> > > values_by_instace_and_time;
        std::vector< std::vector<double> > values_at_time;
        std::vector<double> current_values;

        Eigenvectors(int length_of_D_ij, int num_points_recorded, int num_of_MF_instances);
    };

    int Recieve_D_ij_from_processes(Eigenvectors &eigen, int numProcs, int length_Da_ij, MPI_Status &stat);
    int Calculate_eigenvalues(Eigenvectors &eigen);
    int Calculate_diffusion_coefficients(Eigenvectors &eigen, int numProcs);
    int write_diffusion_coefficients_to_file(Eigenvectors &eigen, std::string filename);

/*******************************************************************/

std::ostream& operator<<(std::ostream& os, const std::vector<double> &v);
std::ostream& operator<<(std::ostream& os, const std::array<double, 3> &v);

#endif