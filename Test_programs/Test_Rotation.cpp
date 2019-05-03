#include "Bfield/class_bfield.h"
#include "Particle/class_particle.h"
#include "Trajectory/class_trajectory.h"
#include "Functions/functions.h"


#define _USE_MATH_DEFINES

constexpr int N_TEST_PARTICLES      =   100;
constexpr int N_RANDOM_MODES        =   500;
constexpr int T_RUN_FOR_YEARS       =   1e5;
constexpr double B_REGULAR_COMP     =   1.0;                                         //microGauss
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
constexpr bool GEN_TURB_GLOB        =   false;


int initialize_init(Initializer &init, int procID);


int main(){


        Initializer init;
    initialize_init(init, 1);

    Ran rng(init.seed);

    Bfield bfield(init, rng);
    Particle particle(init);
    Trajectory trajectory(init);

    std::vector<double> axis, b_hat_prev;
 
    double theta;

    trajectory.initialize_new_GC(particle, bfield, rng);

    b_hat_prev = { 0, 1, 1 };
    bfield.B_hat = {0, 1, -1 };

    //calculate particle.v
        axis = GCT::calculate_vector_rotation(b_hat_prev, bfield.B_hat, theta);     // Axis now holds the normal vector for the rotation

    std::cout << "B_before: " << b_hat_prev << std::endl;
    std::cout << "B_after: " << bfield.B_hat << std::endl << std::endl;

    std::cout << "axis: " << axis << std::endl;
    std::cout << "Theta: " << theta << std::endl << std::endl;
    std::cout << "hat_1: " << trajectory.hat_1 << std::endl;
    std::cout << "a_hat: " << trajectory.a_hat << std::endl;
    

        GCT::rotate_vector_in_plane(trajectory.hat_1, axis, theta);                 // \hat{1} is now rotated to the new plane
        trajectory.a_hat = trajectory.hat_1;                                        // set phase to zero
        GCT::rotate_vector_in_plane(trajectory.a_hat, bfield.B_hat, M_PI);          // rotate by the gyrophase

    std::cout << "\nAfter rotation:\nhat_1: " << trajectory.hat_1 << std::endl;
    std::cout << "a_hat: " << trajectory.a_hat << std::endl;

    std::vector<double> v_perp_hat = GCT::vector_cross_product(trajectory.a_hat, bfield.B_hat);
        GCT::normalize_vector(v_perp_hat);
        particle.v = {                                                              // Add velocity in perp and parallell directions
            trajectory.v_perp*v_perp_hat[0] + trajectory.v_parallell*bfield.B_hat[0],
            trajectory.v_perp*v_perp_hat[1] + trajectory.v_parallell*bfield.B_hat[1],
            trajectory.v_perp*v_perp_hat[2] + trajectory.v_parallell*bfield.B_hat[2],
        };

    return 0;

    
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

  return 0;
}