#include "Functions/functions.h"

#include <boost/mpi.hpp>

/***********
 * These functions are all tested and working for all cases.
***********/

Eigenvectors::Eigenvectors(int length_of_D_ij, int num_points_recorded, int num_of_MF_instances){
    all_Da_ij.resize(num_of_MF_instances, std::vector<double>(num_points_recorded, 0));
    Da_ij.resize(length_of_D_ij);
    D_average.resize(num_points_recorded);

    values_by_instace_and_time.resize(num_of_MF_instances, std::vector< std::vector<double> >(num_points_recorded, std::vector<double>(3, 0)));
    values_at_time.resize(num_points_recorded, std::vector<double>(3,0));
    current_values.resize(3,0);
}

int Recieve_D_ij_from_processes(Eigenvectors &eigen, int numProcs, int length_Da_ij, MPI_Status &stat){
    
    for(int i = 1; i < numProcs; i++){                                  //Start from process 1, as this is process 0

      MPI_Recv(eigen.Da_ij.data(), length_Da_ij, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
      
      eigen.all_Da_ij[stat.MPI_SOURCE - 1] = eigen.Da_ij;
      //D_ij[a] will now hold the matrix for D^(a)_ij

      std::cout << "I recieved D_ij from process " << stat.MPI_SOURCE << std::endl;
    }//End for i

    return 0;
}

int Calculate_eigenvalues(Eigenvectors &eigen){
    
    for(int instance = 0; instance < eigen.all_Da_ij.size(); instance++){
      for(int i = 0; i < eigen.all_Da_ij[instance].size() - 1; i += 7){

        //print_Dij(D_ij[instance], i);
        GCT::calculate_eigenvalues_3x3_sym(eigen.all_Da_ij[instance], i, eigen.current_values);

        
        eigen.values_by_instace_and_time[instance][i/7] = eigen.current_values; 
        

        // Creates the sum for the average diffusion coefficient, over all instances
        eigen.D_average[i/7] += eigen.all_Da_ij[instance][i+6];


      }//End for i               
    }//End for instance
    
    return 0;
}

int Calculate_diffusion_coefficients(Eigenvectors &eigen, int numProcs){
    
    for(int instance = 0; instance < eigen.values_by_instace_and_time.size(); instance++){
      for(int year = 0; year < eigen.values_by_instace_and_time[instance].size(); year++){
        
        for(int d_i = 0; d_i < 3; d_i++){                  //d_i is the eigenvalue d_1, d_2 or d_3 of the matrix D_ij[instance]
          eigen.values_at_time[year][d_i] += eigen.values_by_instace_and_time[instance][year][d_i];
        }
      
      }
    }//Eigenvalues_by_year now holds the sum in the calculation of yearly average eigenvalue


    for(int year = 0; year < eigen.values_at_time.size(); year++){         //Divide sum by number of instances
      for(int d_i = 0; d_i < 3; d_i++){
        eigen.values_at_time[year][d_i] = eigen.values_at_time[year][d_i] / (numProcs - 1);
      }

      eigen.D_average[year] = eigen.D_average[year] / (numProcs - 1);

    }//eigenvalue_by_year now holds the average eigenvalues for each logarithmic year.

    return 0;
}

int write_diffusion_coefficients_to_file(Eigenvectors &eigen, std::string filename){
    
    std::ofstream file;

    GCT::create_directory_to_file(filename);

    file.open(filename);

    if(!file){
        std::cout << "Couldn't open data/eigenvalues.dat. \n";
        throw("Couldn't open eigenvalues.dat");
    }

    for(int year = 0; year < eigen.values_at_time.size(); year++){
        file << eigen.values_at_time[year][0] << ' ' << eigen.values_at_time[year][1] << ' ' << eigen.values_at_time[year][2];
        file << ' ' << eigen.D_average[year] << '\n';  
    }

    std::cout << "Eigenvalues has been written to data/eigenvalues.dat. \n";
    file.close();

    return 0;
}