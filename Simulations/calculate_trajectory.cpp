#include "Trajectory/class_trajectory.h"
#include "Bfield/class_bfield.h" 
#include "Functions/functions.h"

#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>

constexpr int N_TEST_PARTICLES      =   100;
constexpr int N_RANDOM_MODES        =   100;
constexpr int T_RUN_FOR_YEARS       =   1;
constexpr double B_REGULAR_COMP     =   10.0;                                         //microGauss
constexpr double B_TURBULENT_COMP   =   (1/1.0)*B_REGULAR_COMP;                      //microGauss
constexpr double E_TOTAL            =   1e05;                                        //eV
constexpr double LAMBDA_MAX         =   150.0;                                        //pc
constexpr double LAMBDA_MIN         =   0.027;                                       //pc    (0.27 = Rl/10 for B = 4, E = e16)
constexpr double Q_CHARGE           =   1;                                           //# electron charges
constexpr double M_MASS             =   938.2720813;                                 //MeV/c^2
constexpr double ERROR_MAX          =   1.0e-06;                                     
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

    particle.initialize_new_particle(rng);

    std::cout << "vx0: " << particle.v[0] << " vy0: " << particle.v[1] << " vz0: " << particle.v[2] << std::endl;
    std::cout << "pos: " << particle.pos << '\n'; 
    std::cout << "Gamma: " << particle.gamma_l_NR << std::endl;

    //trajectory.Propagate_particle_wtf(particle, init, rng, bfield);
    trajectory.Propagate_particle_wtf(particle, bfield);


    clock_t end = clock();
    double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Program ended in " << elapsed_secs << " seconds \n";
    return 0;
}

int initialize_init(Initializer &init, int procID){
  init.procID = procID;

  std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
  init.seed = 20; //rand() + 1000 * procID;                                 //Sets the seed for the RNG. Set to const 
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



/***********************************************************************************************/
/***********************************************************************************************/


    //trajectory.generate_bfield(trajectory.t, trajectory.Bx, trajectory.By, trajectory.Bz);
    //trajectory.write_turbulence_to_files();
    //trajectory.write_B_to_file();


    //trajectory2.RK4t_step_size_control();
    //trajectory2.RK_BS_wtf(init);


    /*
    for(int i = 0; i < trajectory2.D_ij.size(); i+=7){
        print_Dij(trajectory2.D_ij, i);
    }
*/

 /*       
    trajectory2.RK4t_step_size_control_nw();

    for(int i = 0; i < trajectory2.D_ij.size(); i += 6){
      for(int j = 0; j < 6; j++){

        //i will represent the time in logarithmic years (i.e t_y = 10^(i/6))
        trajectory2.D_ij[i + j] = trajectory2.D_ij[i + j] / (2.0 * init.N * trajectory2.D_ij_time[i/6]);
        //D_ij has now been calculated for one instance of the magnetic field.

      }
    }//End for i

    std::vector<double> eigenvalues_current(3,0);

    for(int i = 0; i < trajectory2.D_ij.size() - 1; i += 6){
           
        calculate_eigenvalues_3x3_sym(trajectory2.D_ij, i, eigenvalues_current);
        cout << eigenvalues_current[0] << ' ' << eigenvalues_current[1] << ' ' << eigenvalues_current[2] << '\n';
    
    }//End for i 
        
 */   
/*
    ofstream test_rng;
    test_rng.open("data/test_rng.dat");

    if(!test_rng) std::cout << "couldn't open test_rng.dat. \n";

    for(int i = 0; i < 10000; i++){
        test_rng << rng.doub() << ' ' << rng.doub() << '\n';
    }

    test_rng.close();
*/

/* 
    ofstream file4;
    file4.open("data/RNG.dat");

    for(int i = 0; i < 10000; i++){
        file4 << trajectory2.ran.doub() << '\n';
    }

    file4.close();
*/