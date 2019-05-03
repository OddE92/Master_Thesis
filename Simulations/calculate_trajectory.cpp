#include "Trajectory/class_trajectory.h"
#include "Bfield/class_bfield.h" 
#include "Functions/functions.h"

#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>
/*
void print_Dij(vector<double> &inVector, int startMatrix){

  cout << setprecision(5);
  cout << endl;

  cout << inVector[startMatrix + 0] << ' ' << inVector[startMatrix + 1] << ' ' << inVector[startMatrix + 2] << endl;
  cout << inVector[startMatrix + 1] << ' ' << inVector[startMatrix + 3] << ' ' << inVector[startMatrix + 4] << endl;
  cout << inVector[startMatrix + 2] << ' ' << inVector[startMatrix + 4] << ' ' << inVector[startMatrix + 5] << endl;

  cout << endl;

  cout << "D_ij average: " << inVector[startMatrix +6] << endl;
  
  cout << endl;

}*/

int main(void){
    clock_t begin = clock();

    Initializer init;
    
    std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
    int seed = 0; //rand() + 1000;

    Ran rng(seed);

    init.procID = 0;

    std::srand(time(0));
    init.seed = rand();                                                             //Sets the seed for the RNG. Set to 0 
                                                                                    //to generate equal results

    init.n_k = 50;                                                                 //#modes used to generate TMF
    
    init.t_end_y = 1e5;                                                             //Set time in years
    init.num_points_recorded = log10(init.t_end_y) * 9 + 1;
    init.t_start = 0.0; 
    init.t_end = 31557600.0 * init.t_end_y; 
    //init.t_end = 1000.0;
    init.dt = 0.1;                                                                  //Set time start and end here, and the timestep

/*
    Set B_0 as the regular magnetic field strength, given in microgauss.
    Set B_rms_turb as the turbulent magnetic field RMS-strength, given in microgauss.
*/
    init.B_0 = 1.0; init.B_rms_turb = 4.0;                              //B_0 is the regular field, B_rms_turb is the RMS of the turb field
    init.lambda_max = 10; init.lambda_min = 0.0027;                     //set these in parsecs
/*
    Set initial conditions. E is the total energy given in eV.
    x0, y0 and z0 are the initial coordinates, given in parsec.
    vx0, vy0 and vz0 are the initial velocities, given as a percentage (must add to 1).
    The speed in each direction is calculated on the basis of the total energy and the
    given percentages.
*/
    init.E = 1.0 * 1e17;                                                
    init.start_pos = { 0, 0, 0 };                        //Initial conditions
    init.start_vel = { 0.9, 0, 1-0.9*0.9};
    init.N = 1;
/* 
    Initialize charge Q, where Q = q*e, q being an integer and e the elementary charge
    Mass is given in MeV/c^2
    for(int i = 0; i < trajectory2.D_ij.size(); i+=7){
        print_Dij(trajectory2.D_ij, i);
    }
*/
    init.q = 1.0; init.m = 938.2720813;                                 //Set the charge and mass
                                                                        
    init.max_err = 1.0e-18;                                             //Set min and max error
    init.min_err = 1.0e-08;

    init.generate_turbulence = false;                                    //This does decide if you generate a turbulence or not
    init.GCT = false;
    
    //Trajectory trajectory1(init);
    Bfield bfield(init, rng);
    Particle particle(init);
    Trajectory trajectory(init);

    particle.initialize_new_particle(rng);
    particle.v = { 0.9*c, 0, std::sqrt(1-0.9*0.9)*c};
    std::cout << "vx0: " << particle.v[0] << " vy0: " << particle.v[1] << " vz0: " << particle.v[2] << std::endl;   

    trajectory.Propagate_particle_wtf(bfield, particle);

    //trajectory.generate_bfield(trajectory.t, trajectory.Bx, trajectory.By, trajectory.Bz);
    //trajectory.write_turbulence_to_files();
    //trajectory.write_B_to_file();


    //trajectory2.RK4t_step_size_control();
    //trajectory2.RK_BS_wtf(init);

    trajectory.R_Larmor = 1.0810076e-15 * (particle.v[0] / c) * (particle.E / (bfield.B_0 + bfield.B_rms_turb*0)); 

    std::cout << "R_l: " << trajectory.R_Larmor << std::endl;
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
    clock_t end = clock();
    double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Program ended in " << elapsed_secs << " seconds \n";

}