#include "Bfield/class_bfield.h"
#include "NR3/ran.h"
#include "Particle/class_particle.h"
#include "Trajectory/class_trajectory.h"
#include "Functions/functions.h"
#include "Initializer/initializer.h"

#include <boost/mpi.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

constexpr int N_TEST_PARTICLES      =   100;
constexpr int N_RANDOM_MODES        =   50;
constexpr int T_RUN_FOR_YEARS       =   1e5;
constexpr double B_REGULAR_COMP     =   0.0;                                         //microGauss
constexpr double B_TURBULENT_COMP   =   4.0;                                         //microGauss
constexpr double E_TOTAL            =   1e18;                                        //eV
constexpr double LAMBDA_MAX         =   150.0;                                        //pc
constexpr double LAMBDA_MIN         =   0.027;                                       //pc    (0.27 = Rl/10 for B = 4, E = e16)
constexpr double Q_CHARGE           =   1;                                           //# electron charges
constexpr double M_MASS             =   938.2720813;                                 //MeV/c^2
constexpr double ERROR_MAX          =   1.0e-05;                                     
constexpr double ERROR_MIN          =   1.0e-08;
const int NUM_POINTS_RECORDED       =   log10(T_RUN_FOR_YEARS) * 9 + 1;
const int D_IJ_LENGTH               =   NUM_POINTS_RECORDED * 7;
constexpr bool GEN_TURB_GLOB        =   true;


int initialize_init(Initializer &init, int procID);


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int procID;
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);
  //procID = 1;       //For testing

 /******************************************/ 
 /************* MAIN PROCESSES *************/
  if(procID == 0){                                                      // Process to finalize calculations
    clock_t begin = clock();

    MPI_Status stat;
    
    Eigenvectors eigen(D_IJ_LENGTH, NUM_POINTS_RECORDED, numProcs - 1);


    Recieve_D_ij_from_processes(eigen, numProcs, D_IJ_LENGTH, stat);
    // eigen.all_Da_ij now holds all the symmetric matrices, sorted by instance 
    // i.e. eigen.all_Da_ij[0] holds the matrix for the first MF-instance


    Calculate_eigenvalues(eigen);
    // eigen.values_by_instace_and_time now holds the eigenvalues
    // The sum of average diffusion coefficients have also been calculated


    Calculate_diffusion_coefficients(eigen, numProcs);
    // Sums the eigenvalues over all instances and divides by number of instaces
    // to calculate the diffusion coefficients d_1, d_2 and d_3


    std::string filename = GCT::generate_unique_filename_eigenvalues(E_TOTAL, LAMBDA_MAX, procID);

    write_diffusion_coefficients_to_file(eigen, filename);


    //This should now have calculated the eigenvalues.


    clock_t end = clock();
    double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Rank 0 ended in " << elapsed_secs << " seconds \n";
    std::cout << "Values are: E = " << E_TOTAL << ", L_max = " << LAMBDA_MAX << std::endl;



 /*************************************************/    
 /************* CALCULATION PROCESSES *************/
  }else{            


    double ranPhi, ranTheta, ranTot;
    double percentCounter = 0.099999;

    Initializer init;

    initialize_init(init, procID);

    Ran rng(15321 + 100*procID + init.seed);                            //To generate directions
    Bfield bfield(init, rng);
    Particle particle(init);
    Trajectory trajectory(init);
    
    std::ofstream file;                                                 //To hold position of particles
    std::string filename = GCT::generate_unique_filename_positions(bfield, particle, procID);
    GCT::create_directory_to_file(filename);
    file.open(filename);
    
    if(!file){
      std::cout << "Couldn't open " << filename << std::endl;
      throw("Couldn't open " + filename);
    }
    

    for(int i = 0; i < init.N; i++){                                    //i < number of particles to test

      std::cout << "When do I get here?\n";
      particle.initialize_new_particle(rng);
      std::cout << "What about here?\n";

      trajectory.Propagate_particle(particle, bfield);

      std::cout << "I did a particle\n";

      //Print progress in %
      if( static_cast<double>(i+1)/init.N - percentCounter > __DBL_EPSILON__ ){
        std::cout << "Progress rank " << procID << ": " << (static_cast<double>(i+1)/init.N) * 100 << "\% \n";
        percentCounter += 0.1;
      }


      trajectory.write_positions_to_file(file);

     //After this loop, D_ij holds the sum in the equation for D_ij
    }//End for particle i


    trajectory.calculate_D_ij(init.N);
    
    
    std::cout << "Number of coordinate sets in r_vect: " << init.N * trajectory.r_vect.size() / 3.0 << std::endl;
    MPI_Send(trajectory.D_ij.data(), trajectory.D_ij.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  
  }
 /************* END CALCULATION PROCESSES *************/

  MPI_Finalize();
  return 0;
}

/************* END GENERATE SAMPLES *************/
/************************************************/



int initialize_init(Initializer &init, int procID){
  init.procID = procID;

  std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
  init.seed = 1000 * procID;                                 //Sets the seed for the RNG. Set to const 
                                                                      //to generate equal results

  init.n_k = N_RANDOM_MODES;                                          //#modes used to generate TMF
  init.N = N_TEST_PARTICLES;                                          //Number of test-particles
    
  init.t_end_y = T_RUN_FOR_YEARS;                                     //Set time in years
  init.t_start = 0.0;                                                 //Set time start and end here, and the timestep
  init.t_end = 31557600.0 * init.t_end_y; 
  init.dt = 0.1;        

  init.B_0 = B_REGULAR_COMP; init.B_rms_turb = B_TURBULENT_COMP;      //B_0 is the regular field, B_rms_turb is the RMS of the turb field

  init.E = E_TOTAL;                                                   //Energy in eV

  init.lambda_max = LAMBDA_MAX; init.lambda_min = LAMBDA_MIN;         //Wavelength in pc

  init.q = Q_CHARGE; init.m = M_MASS;                                 //Set the charge and mass
                                                                        
  init.max_err = ERROR_MAX;                                           //Set min and max error
  init.min_err = ERROR_MIN;
                              
  init.generate_turbulence = GEN_TURB_GLOB;                           //This decides if you generate a turbulence or not
    
  init.num_points_recorded = NUM_POINTS_RECORDED;                     //Sets length of diffusion tensor-vector  

  init.start_pos = { 0, 0, 0 };                                       //Set here as default. Should be randomly chosen in the program
  init.start_vel = { 1, 0, 0 };

  init.GCT = false;
  return 0;
}

