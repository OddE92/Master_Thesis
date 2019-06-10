#ifndef CLASS_TRAJECTORY
#define CLASS_TRAJECTORY

#include "NR3/nr3.h"
#include "Initializer/initializer.h"
#include "Particle/class_particle.h"
#include "Bfield/class_bfield.h"
#include "Guiding_center/class_guiding_center.h"
#include "Functions/functions.h"
#include "Constants/constants.h"

#include <vector>
#include <array>

class Trajectory{
    public:
      // For holding diffusion tensor and particle positions
        std::vector<double> D_ij, D_ij_time;
        std::vector<double> r_vect;
        
        double t, t_start, t_end, dt;
        double min_err, max_err;
        double unit_coeff;
        double R_Larmor;
        double timestep;

       // RK-solver variables
        int nvar;
        VecDoub ystart;

       // Functions  
        //int Propagate_particle(Particle &particle, Initializer &init, Ran &rng);
        int Propagate_particle(Particle &particle, Bfield &bfield);
        //int Propagate_particle_wtf(Particle &particle, Initializer &init, Ran &rng, Bfield &bfield);
        int Propagate_particle_wtf(Particle &particle, Bfield &bfield);
        int Propagate_GC(Bfield &bfield, Particle &particle, Guiding_Center &GC, Initializer &init, Ran &rng);
        int Propagate_GC_wtf(Bfield &bfield, Particle &particle, Guiding_Center &GC, Initializer &init, Ran &rng);
        
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

    std::array<double, 3> B, r;

    Bfield bfield;

    Rhs_lorentz_equation(Particle &particle){ 
        unit_coeff = 8.987 * particle.q / (particle.gamma_l * particle.m); 
    }
 
    Rhs_lorentz_equation(Particle &particle, Bfield &bbfield) : bfield(bbfield){  
        unit_coeff = 8.987 * particle.q / (particle.gamma_l * particle.m);      //Coefficient for the units when dv/dt = unit_coeff * V x B
    }
 
    /*
        y[0] = x;  y[1] = y;  y[2] = z;
        y[3] = vx; y[4] = vy; y[5] = vz;
        y[6] = Bx; y[7] = By; y[8] = Bz;
    */

    void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydx){

            // Generates the bfield at the point, then calculates the new velocity and positions

            r = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };

            bfield.generate_bfield_at_point(t, B, r);            
     
        dydx[0] = y[3];                         // dydx[0] = vx
        dydx[1] = y[4];                         // dydx[1] = vy
        dydx[2] = y[5];                         // dydx[2] = vz
     
        dydx[3] = unit_coeff * (y[4] * B[2] - y[5] * B[1]);    // ax = vy*Bz - vz*By
        dydx[4] = unit_coeff * (y[5] * B[0] - y[3] * B[2]);    // ay = vz*Bx - vx*Bz
        dydx[5] = unit_coeff * (y[3] * B[1] - y[4] * B[0]);    // az = vx*By - vy*Bx
     /*
        dydx[3] = unit_coeff * (y[4] * y[8] - y[5] * y[7]);    // ax = vy*Bz - vz*By
        dydx[4] = unit_coeff * (y[5] * y[6] - y[3] * y[8]);    // ay = vz*Bx - vx*Bz
        dydx[5] = unit_coeff * (y[3] * y[7] - y[4] * y[6]);    // az = vx*By - vy*Bx

        dydx[6] = 0;
        dydx[7] = 0;
        dydx[8] = 0;
     */
    } 

};

struct GC_equation{

    /*****
     * 
     *  The guiding center equation is more complex than the Lorentz-equations, and 
     *  requires more info at each step. Thus it has its own particle and GC 
     *  as well as a bfield. It also keeps some extra doubles and arrays to
     *  hold all values needed in the calculations.
     * 
     *  operator() is ound in class_trajectory.cpp
     * 
     *  The equations are taken from:
     *      Cary, John R. and Brizard, Alain J.
     *      Hamiltonian Theory of Guiding-Center Motion
     *      https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.693
     * 
     *  Equations:
     *      dy/dt[0-2]  = (3.13)
     *      dy/dt[3]    = (3.12)
     *      dy/dt[4]    = (3.7)
     * 
    *****/

 /*
    y[0-2]  =   GC_position
    y[3]    =   u
    y[4]    =   gyrophase
 */

    Bfield bfield;
    Particle particle;
    Guiding_Center GC;

    double q, m, theta;

    std::array<double, 3> r, E_eff_cross_B_hat, b_hat_prev, axis;

    int calculate_particle_v();

    double B_eff_dot_E_eff, B_eff_parallell, B_eff_amp;

    GC_equation(Particle &pparticle, Bfield &bbfield, Guiding_Center &GGC) : particle(pparticle), bfield(bbfield), GC(GGC){
        q = particle.q; m = particle.m;
        b_hat_prev = bbfield.B_hat;
    }

    void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydx);

};




#endif  