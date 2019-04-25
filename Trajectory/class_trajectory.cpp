#include "Trajectory/class_trajectory.h"

#include "Integrator/odeint.h"
#include "Integrator/stepperbase.h"
#include "Integrator/stepperBS.h"
#include "Integrator/stepperdopr5.h"
#include "Integrator/stepperdopr853.h"

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





int Trajectory::Propagate_particle(Bfield &bfield, Particle &particle){
    double hmin = 0.0;                  // Needed in Odeint
    this->t = t_start;                  // Reset t before propagating the particle

    R_Larmor = 1.0810076e-15 * (particle.v_total / c) * (particle.E / (bfield.B_0 + bfield.B_rms_turb)); 

    ystart[0] = particle.pos[0]; ystart[1] = particle.pos[1]; ystart[2] = particle.pos[2];
    ystart[3] = particle.v[0]; ystart[4] = particle.v[1]; ystart[5] = particle.v[2];

    Output out;

    Rhs_lorentz_equation r(particle);

    Odeint<StepperBS<Rhs_lorentz_equation> > ode(ystart, t_start, t_end, max_err, max_err, dt, hmin, out, r);

    ode.integrate(this->D_ij, this->D_ij_time, this->r_vect, bfield, particle, this->t_end);

    return 0;
}

int Trajectory::Propagate_particle_wtf(Bfield &bfield, Particle &particle){         // Used to save a single particle trajectory
    double hmin = 0.0;

    R_Larmor = 1.0810076e-15 * (particle.v_total / c) * (particle.E / (bfield.B_0 + bfield.B_rms_turb)); 

    ystart[0] = particle.pos[0]; ystart[1] = particle.pos[1]; ystart[2] = particle.pos[2];
    ystart[3] = particle.v[0]; ystart[4] = particle.v[1]; ystart[5] = particle.v[2];

    Output out;

    Rhs_lorentz_equation r(particle);

    Odeint<StepperBS<Rhs_lorentz_equation> > ode(ystart, t_start, t_end, max_err, max_err, dt, hmin, out, r);
    
    ode.integrate(this->D_ij, this->D_ij_time, this->r_vect, bfield, particle, this->t_end);

    std::cout << "NOK: " << ode.nok << "; Nbad: " << ode.nbad << std::endl;

    std::ofstream file;

    file.open("Data/RK853.dat");
    std::cout << "Ant. punkter:" << out.count << std::endl;
    for(int i = 0; i < out.count; i++){
        file << out.ysave[0][i]*mtopc << ' ' << out.ysave[1][i]*mtopc << ' ' << out.ysave[2][i]*mtopc << '\n';
    }

    file.close();

    return 0;
}




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

    ystart.resize(nvar);
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