#include "Trajectory/class_trajectory.h"

#include "Integrator/odeint.h"
#include "Integrator/stepperbase.h"
#include "Integrator/stepperBS.h"
#include "Integrator/stepperdopr5.h"
#include "Integrator/stepperdopr853.h"
#include "Functions/functions.h"

#include <iostream>

int Trajectory::calculate_D_ij(int num_particles){

    for(int i = 0; i < this->D_ij.size(); i += 7){
      for(int j = 0; j < 7; j++){

        this->D_ij[i + j] = this->D_ij[i + j] / (2.0 * num_particles * this->D_ij_time[i/7]);      //std::pow(10, i/6)
        //D_ij has now been calculated for one instance of the magnetic field.

        if(j == 6) this->D_ij[i + j] = this->D_ij[i + j] / 3;     //To calculate the average diffusion coefficient

      }
    }//End for i
    return 0;
}



int Trajectory::write_positions_to_file(std::ofstream &file){
    //Write position of particles to file
    for(int j = 0; j < this->r_vect.size(); j += 3){
      file << this->r_vect[j+0] << ' ' << this->r_vect[j+1] << ' ' << this->r_vect[j+2] << ' ';
      file << pow(this->r_vect[j+0], 2) + pow(this->r_vect[j+1], 2) + pow(this->r_vect[j+2], 2) << '\n';
    }
    return 0;
}



    // This function sets the given initial values, then propagates the passed particle
    // The bfield, particle and GC is copied into GC_equation r
int Trajectory::Propagate_GC(Bfield &bfield, Particle &particle, Guiding_Center &GC, Initializer &init, Ran &rng){
    double hmin = 0.0;                  // Needed in Odeint
    this->t = t_start;                  // Reset t before propagating the particle

    Output out;                         // default constructor disregards all results (i.e no saves)
  
  // Set initial values 
    ystart[0] = GC.GC_position[0]; ystart[1] = GC.GC_position[1]; ystart[2] = GC.GC_position[2];
    ystart[3] = GC.u;
    ystart[4] = GC.gyrophase;

    GC_equation r(particle, bfield, GC);

    Odeint<StepperBS<GC_equation> > ode(ystart, t_start, t_end, max_err, max_err, dt, hmin, out, r);

    // r can be accessed inside integrate_GC as "derivs", i.e "derivs.GC.timestep" gives r.GC.timestep.

    ode.integrate_GC(*this, this->t_end);

    return 0;
}

    // Same function as above, but writes the result to a file
int Trajectory::Propagate_GC_wtf(Bfield &bfield, Particle &particle, Guiding_Center &GC, Initializer &init, Ran &rng){
    double hmin = 0.0;                  // Needed in Odeint
    this->t = t_start;                  // Reset t before propagating the particle  

    Output out(-1);                     // out(-anyNumber) saves the values of y at each point

    ystart[0] = GC.GC_position[0]; ystart[1] = GC.GC_position[1]; ystart[2] = GC.GC_position[2];
    ystart[3] = GC.u;
    ystart[4] = GC.gyrophase;

    GC_equation r(particle, bfield, GC);

    Odeint<StepperBS<GC_equation> > ode(ystart, t_start, t_end, max_err, max_err, dt, hmin, out, r);
    
    ode.integrate_GC(*this, this->t_end);

    std::ofstream file;

    file.open("Data/trajectory_GC.dat");
    std::cout << "Ant. punkter:" << out.count << std::endl;
    for(int i = 0; i < out.count; i++){
        file << out.ysave[0][i]*GCT::mtopc << ' ' << out.ysave[1][i]*GCT::mtopc << ' ' << out.ysave[2][i]*GCT::mtopc << '\n';
    }

    file.close();


    return 0;
}





    // Same as above, but for the standard Lorentz force.
int Trajectory::Propagate_particle(Particle &particle, Initializer &init, Ran &rng){
    double hmin = 0.0;                  // Needed in Odeint
    this->t = t_start;                  // Reset t before propagating the particle

    ystart[0] = particle.pos[0]; ystart[1] = particle.pos[1]; ystart[2] = particle.pos[2];
    ystart[3] = particle.v[0]; ystart[4] = particle.v[1]; ystart[5] = particle.v[2];

    Output out;

    Rhs_lorentz_equation r(particle, init, rng);

    Odeint<StepperBS<Rhs_lorentz_equation> > ode(ystart, t_start, t_end, max_err, max_err, dt, hmin, out, r);

    ode.integrate(this->D_ij, this->D_ij_time, this->r_vect, particle, this->t_end);

    return 0;
}
    // Same as above, writes to file
int Trajectory::Propagate_particle_wtf(Particle &particle, Initializer &init, Ran &rng){         // Used to save a single particle trajectory
    double hmin = 0.0;

    ystart[0] = particle.pos[0]; ystart[1] = particle.pos[1]; ystart[2] = particle.pos[2];
    ystart[3] = particle.v[0]; ystart[4] = particle.v[1]; ystart[5] = particle.v[2];

    Output out(-1);

    Rhs_lorentz_equation r(particle, init, rng);

    Odeint<StepperBS<Rhs_lorentz_equation> > ode(ystart, t_start, t_end, max_err, max_err, dt, hmin, out, r);
    
    ode.integrate(this->D_ij, this->D_ij_time, this->r_vect, particle, this->t_end);

    std::cout << "NOK: " << ode.nok << "; Nbad: " << ode.nbad << std::endl;

    std::ofstream file;

    file.open("Data/trajectory.dat");
    std::cout << "Ant. punkter:" << out.count << std::endl;
    for(int i = 0; i < out.count; i++){
        file << out.ysave[0][i]*GCT::mtopc << ' ' << out.ysave[1][i]*GCT::mtopc << ' ' << out.ysave[2][i]*GCT::mtopc << '\n';
    }

    file.close();

    return 0;
}


    // Initializers
Trajectory::Trajectory(Initializer &init){
    
    t = 0.0;
    t_start = init.t_start;
    t_end = init.t_end;
    dt = init.dt;

    max_err = init.max_err;
    min_err = init.min_err;

    D_ij.resize(init.num_points_recorded * 7);            // 6 parts for the matrix and one for r^2 (D_average)
    D_ij_time.resize(init.num_points_recorded);

    r_vect.resize(init.num_points_recorded * 3);          // Holds position at each recorded point

    if(init.GCT){ nvar = 5;
    }else {nvar = 6;}

    ystart.resize(nvar);                                  // Used to set initial values
}

Trajectory::Trajectory(){

    std::cout << "Trajectory: Default constructar called. \n\n";

    t = 0.0;
    t_start = 0.0;
    t_end = 1.0e5 * 31557600.0;                 // 10^5 years in [s]
    dt = 0.1;

    max_err = 1.0e-05;
    min_err = 1.0e-08;

   // Default is 1e5 years
    D_ij.resize(46 * 7); 
    D_ij_time.resize(46);

    r_vect.resize(46 * 3);
}



/**********
 * 
 *  Operator() for the guiding center equation.
 *  Currently under work
 * 
**********/


void GC_equation::operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydx){

        r = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
    std::cout << "r: " << r << std::endl;
        bfield.generate_bfield_at_point(t, r);
        
        GC.GC_velocity = { y[3]*bfield.B_hat[0], y[3]*bfield.B_hat[1], y[3]*bfield.B_hat[2] };
   
  /********** CALCULATE PARTICLE.V **********/
    std::cout << "B_hat_prev: " << b_hat_prev << std::endl;
    std::cout << "B_hat: " << bfield.B_hat << std::endl;
        axis = GCT::calculate_vector_rotation(b_hat_prev, bfield.B_hat, theta);     // Axis now holds the normal vector for the rotation
    std::cout << "Axis, prenorm: " << axis << std::endl;
        GCT::normalize_vector(axis);
    std::cout << "theta: " << theta << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    
        GCT::rotate_vector_in_plane(GC.hat_1, axis, theta);                 // \hat{1} is now rotated to the new plane
        GC.a_hat = GC.hat_1;                                                // set phase to zero
    std::cout << "hat_1: " << GC.hat_1 << std::endl;
        GCT::rotate_vector_in_plane(GC.a_hat, bfield.B_hat, GC.gyrophase);          // rotate by the gyrophase
    std::cout << "B_hat: " << bfield.B_hat << std::endl;
    std::cout << "a_hat: " << GC.a_hat << std::endl;
        GC.v_perp_hat = GCT::vector_cross_product(GC.a_hat, bfield.B_hat);
    std::cout << "v_perp_hat: " << GC.v_perp_hat << std::endl;
        GCT::normalize_vector(GC.v_perp_hat);
        particle.v = {                                                                  // Add velocity in perp and parallell directions
            GC.v_perp*GC.v_perp_hat[0] + GC.v_parallell*bfield.B_hat[0],
            GC.v_perp*GC.v_perp_hat[1] + GC.v_parallell*bfield.B_hat[1],
            GC.v_perp*GC.v_perp_hat[2] + GC.v_parallell*bfield.B_hat[2],
        };
  /********** END CALCULATE PARTICLE.V **********/

        std::cout << "GC_pos: " << y[0] << ' ' << y[1] << ' ' << y[2] << std::endl;
        std::cout << "Particl.v: " << particle.v << std::endl;
        std::cout << "gyrophase: " << y[4] << std::endl;
        std::cout << "u: " << y[3] << std::endl;
 
        bfield.calculate_partial_b_hat(bfield.B, GC.GC_velocity, GC.GC_position, t, GC.timestep);

 //std::cout << "B_partial: \n" << bfield.B_hat_partial[0] << std::endl << bfield.B_hat_partial[1] << std::endl << bfield.B_hat_partial[2] << std::endl;


        bfield.calculate_B_effective(GC.GC_velocity, r, y[3]);
        bfield.calculate_E_effective(particle, GC.v_perp);                              // Assume v_perp = const
     std::cout << "B_eff: " << bfield.B_effective << std::endl;
     std::cout << "E_eff: " << bfield.E_effective << std::endl << std::endl;
        E_eff_cross_B_hat = GCT::vector_cross_product(bfield.E_effective, bfield.B_hat);

        B_eff_dot_E_eff = GCT::vector_dot_product(bfield.B_effective, bfield.E_effective);

        B_eff_parallell = GCT::vector_dot_product(bfield.B_hat, bfield.B_effective);

        B_eff_amp = GCT::vector_amplitude(bfield.B_effective);


            if(std::isnan(bfield.B_effective[0])){ 
                std::cout << "B_eff is nan\n";
                throw("Error");
                }

        for(int i = 0; i < 3; i++) dydx[i] =                                                        // dydx[0-2] = GC_v[0-2]
            (y[3]/B_eff_parallell)*bfield.B_effective[i] + (1/B_eff_parallell)*E_eff_cross_B_hat[i];    

        dydx[3] = (1/GCT::mq) * (q / m) * (1 / B_eff_parallell) * B_eff_dot_E_eff;                  // dydx[3] = d/dt u

        dydx[4] = (1/GCT::mq) * (q / (m*GCT::c)) * B_eff_amp;                                       // dydx[4] = d/dt gyrophase

 /*
    y[0-2]  =   GC_position
    y[3]    =   u
    y[4]    =   gyrophase

 */
    }