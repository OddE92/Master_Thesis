#include "Trajectory/class_trajectory.h"
#include "Bfield/class_bfield.h" 
#include "Functions/functions.h"
#include "Particle/class_particle.h"

#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>

constexpr int N_TEST_PARTICLES      =   100;
constexpr int N_RANDOM_MODES        =   500;
constexpr int T_RUN_FOR_YEARS       =   1e5;
constexpr double B_REGULAR_COMP     =   1.0;                                         //microGauss
constexpr double B_TURBULENT_COMP   =   4.0;                                         //microGauss
constexpr double E_TOTAL            =   1e17;                                        //eV
constexpr double LAMBDA_MAX         =   150.0;                                        //pc
constexpr double LAMBDA_MIN         =   0.027;                                       //pc    (0.27 = Rl/10 for B = 4, E = e16)
constexpr double Q_CHARGE           =   1;                                           //# electron charges
constexpr double M_MASS             =   938.2720813;                                 //MeV/c^2
constexpr double ERROR_MAX          =   1.0e-05;                                     
constexpr double ERROR_MIN          =   1.0e-08;
const int NUM_POINTS_RECORDED       =   log10(T_RUN_FOR_YEARS) * 9 + 1;
const int D_IJ_LENGTH               =   NUM_POINTS_RECORDED * 7;
constexpr bool GEN_TURB_GLOB        =   false;

int initialize_init(Initializer &init, int procID);

int main(void){
    clock_t begin = clock();

    Initializer init;

    initialize_init(init, 1);

    Ran rng(init.seed);
    Bfield bfield(init, rng);
    Particle particle(init);
    Trajectory trajectory(init);


    trajectory.initialize_new_GC(particle, bfield, rng);
 std::cout << "GC initial velocity: " << trajectory.GC_velocity << std::endl;


    trajectory.Propagate_GC_wtf(bfield, particle);



    std::cout << "R_l: " << trajectory.R_Larmor << std::endl;

    clock_t end = clock();
    double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Program ended in " << elapsed_secs << " seconds \n";

}

int initialize_init(Initializer &init, int procID){
  init.procID = procID;

  std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
  init.seed = rand() + 1000 * procID;                                 //Sets the seed for the RNG. Set to const 
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

  init.GCT = true;
  return 0;
}