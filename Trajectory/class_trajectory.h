#ifndef CLASS_TRAJECTORY
#define CLASS_TRAJECTORY

#include "NR3/nr3.h"
#include "Initializer/initializer.h"
#include "Particle/class_particle.h"
#include "Bfield/class_bfield.h"

#include <vector>

class Trajectory{
    public:

        std::vector<double> D_ij, D_ij_time;
        std::vector<double> r_vect;
        double t, t_start, t_end, dt;
        double min_err, max_err;
        double unit_coeff;
        double R_Larmor;
        const double mtopc = 3.24078e-17;

       // RK-solver variables
        const int nvar = 9;
        VecDoub ystart;

       // Functions  
        int Propagate_particle(Bfield &bfield, Particle &particle);
        int Propagate_particle_wtf(Bfield &bfield, Particle &particle);
        
        int calculate_D_ij(int num_particles);

        int write_positions_to_file(std::ofstream &file);

        Trajectory();
        Trajectory(Initializer &init);
        ~Trajectory(){};
};

/*************************************************/
/********** Equation of Motion Functors **********/

struct Rhs_lorentz_equation{

    double unit_coeff;

    Rhs_lorentz_equation(Particle &particle){ 
        unit_coeff = 8.987 * particle.q / (particle.gamma_l * particle.m);                 //Coefficient for the units when dv/dt = unit_coeff * V x B
    }

    /*
        y[0] = x;  y[1] = y;  y[2] = z;
        y[3] = vx; y[4] = vy; y[5] = vz;
        y[6] = Bx; y[7] = By; y[8] = Bz;
    */

    void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydx){

        dydx[0] = y[3];                         // dydx[0] = vx
        dydx[1] = y[4];                         // dydx[1] = vy
        dydx[2] = y[5];                         // dydx[2] = vz

        dydx[3] = unit_coeff * (y[4] * y[8] - y[5] * y[7]);    // ax = vy*Bz - vz*By
        dydx[4] = unit_coeff * (y[5] * y[6] - y[3] * y[8]);    // ay = vz*Bx - vx*Bz
        dydx[5] = unit_coeff * (y[3] * y[7] - y[4] * y[6]);    // az = vx*By - vy*Bx

        dydx[6] = 0;                            // 6, 7 and 8 are to keep the B-field
        dydx[7] = 0;                            // it is updated at each step in Odeint::integrate();
        dydx[8] = 0;

    } 

};
/*
struct gc_equation{

    gc_equation(Trajectory_initializer &init){}

    void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydx){

        double B_par = sqrt(y[6]*y[6] + y[7]*y[7] + y[8]*y[8]);

        dydx[0] = y[3];                         // dydx[0] = vx
        dydx[1] = y[4];                         // dydx[1] = vy
        dydx[2] = y[5];                         // dydx[2] = vz

        dydx[3] = u * y[6] / B_par;             // dX/dx = u * B_x / B_par
        dydx[4] = u * y[7] / B_par;             // dX/dy = u * B_y / B_par
        dydx[5] = u * y[8] / B_par;             // dX/dz = u * B_z / B_par

        dydx[6] = 0;                            // 6, 7 and 8 are to keep the B-field
        dydx[7] = 0;                            // it is updated at each step in Odeint::integrate();
        dydx[8] = 0;
    }

};
*/



#endif