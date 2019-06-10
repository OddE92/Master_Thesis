#include "Trajectory/class_trajectory.h"
#include "Bfield/class_bfield.h" 
#include "Functions/functions.h"
#include "Particle/class_particle.h"
#include "Guiding_center/class_guiding_center.h"

#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random>

constexpr int RUN_N_SIMULATIONS     =   1000;

constexpr int N_TEST_PARTICLES      =   100;
constexpr int N_RANDOM_MODES        =   100;
constexpr int T_RUN_FOR_YEARS       =   5e3;
constexpr double B_REGULAR_COMP     =   10.0;                                         //microGauss
constexpr double B_TURBULENT_COMP   =   (1/100.0)*B_REGULAR_COMP;                       //microGauss
constexpr double E_TOTAL            =   1e18;                                        //eV
constexpr double LAMBDA_MAX         =   150.0;                                        //pc
constexpr double LAMBDA_MIN         =   0.027;                                       //pc    (0.27 = Rl/10 for B = 4, E = e16)
constexpr double Q_CHARGE           =   1;                                           //# electron charges
constexpr double M_MASS             =   938.2720813;                                 //MeV/c^2
constexpr double ERROR_MAX          =   1.0e-10;                                     
constexpr double ERROR_MIN          =   1.0e-08;
const int NUM_POINTS_RECORDED       =   log10(T_RUN_FOR_YEARS) * 9 + 1;
const int D_IJ_LENGTH               =   NUM_POINTS_RECORDED * 7;
constexpr bool GEN_TURB_GLOB        =   true;

int initialize_init(Initializer &init, int procID);

int main(int argc, char* argv[]){

  MPI_Init(&argc, &argv);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int procID;
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);

  Initializer init;

  initialize_init(init, procID);

  if(procID == 0){
    //This thread handles the GCT-method

    std::ofstream file;
    int Bpercent = static_cast<int>((B_TURBULENT_COMP / B_REGULAR_COMP)*100 );
    std::string filename = GCT::generate_unique_filename_compare(init.E, Bpercent, procID);
    file.open(filename);

    if(!file) std::cout << "couldn't open the file (thread 0)\n";

    for(int n = 0; n < RUN_N_SIMULATIONS; n++){

      std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

      Ran rng(init.seed);
      Bfield bfield(init, rng);
      Particle particle(init);
      Guiding_Center GC(init);
      Trajectory trajectory(init);
  
      GC.initialize_new_GC(particle, bfield, rng, init.t_start);

      trajectory.Propagate_GC(bfield, particle, GC, init, rng);
  
      std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

      std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t2 - t1);

      file << GC.GC_position << ' ' << GC.R_Larmor << ' ' << time_span.count() << std::endl;

      init.seed++;

      if(n % RUN_N_SIMULATIONS == 0)
      std::cout << "Thread " << procID << ": " << static_cast<int>(100*(n/RUN_N_SIMULATIONS)) << std::endl;

    }//end for

    file.close();

  }
  else if(procID == 1){
    init.GCT = false;

    std::ofstream file;
    int Bpercent = static_cast<int>((B_TURBULENT_COMP / B_REGULAR_COMP)*100 );
    std::string filename = GCT::generate_unique_filename_compare(init.E, Bpercent, procID);
    file.open(filename);

    if(!file) std::cout << "couldn't open the file (thread 1)\n";

    for(int n = 0; n < RUN_N_SIMULATIONS; n++){

      std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

      Ran rng(init.seed);
      Bfield bfield(init, rng);
      Particle particle(init);
      Trajectory trajectory(init);

      particle.initialize_new_particle(rng);

      trajectory.Propagate_particle(particle, bfield);

      std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

      std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t2 - t1);

      file << particle.pos << ' ' << time_span.count() << std::endl;

      init.seed++;

      if(n % RUN_N_SIMULATIONS == 0)
      std::cout << "Thread " << procID << ": " << static_cast<int>(100*(n/RUN_N_SIMULATIONS)) << std::endl;

      

    }//end for

    file.close();

  }
  else{
    std::cout << "Too many threads, smartass. This is thread " << procID << ".\n";
  }

  MPI_Finalize();
  return 0;
}//end main



int initialize_init(Initializer &init, int procID){
  init.procID = procID;

  std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
  init.seed = 0;//rand() + 1000 * procID;                                 //Sets the seed for the RNG. Set to const 
                                                                      //to generate equal results

  init.n_k = N_RANDOM_MODES;                                          //#modes used to generate TMF
  init.N = N_TEST_PARTICLES;                                          //Number of test-particles
    
  init.t_end_y = T_RUN_FOR_YEARS;                                     //Set time in years
  init.t_start = 0.0;                                                 //Set time start and end here, and the timestep
  init.t_end = 31557600.0 * init.t_end_y; 
  init.dt = 10000.0;        

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

  init.GCT = true;
  return 0;
}