#ifndef CLASS_TRAJECTORY
#define CLASS_TRAJECTORY

#include "NR3/nr3.h"
#include "Initializer/initializer.h"
#include "Particle/class_particle.h"
#include "Bfield/class_bfield.h"
#include "Guiding_center/class_guiding_center.h"

#include <vector>

class Trajectory{
    public:

        std::vector<double> D_ij, D_ij_time;
        std::vector<double> r_vect;
        
        double t, t_start, t_end, dt;
        double min_err, max_err;
        double unit_coeff;
        //double R_Larmor;
        const double mtopc = 3.24078e-17;

       // RK-solver variables
        int nvar;
        VecDoub ystart;

       // Functions  
        int Propagate_particle(Bfield &bfield, Particle &particle);
        int Propagate_particle_wtf(Bfield &bfield, Particle &particle);
        int Propagate_GC(Bfield &bfield, Particle &particle, Guiding_Center &GC);
        int Propagate_GC_wtf(Bfield &bfield, Particle &particle, Guiding_Center &GC);
        
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
        unit_coeff = 8.987 * particle.q / (particle.gamma_l * particle.m);      //Coefficient for the units when dv/dt = unit_coeff * V x B
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

struct GC_equation{

    double q, m;
    double qm_coeff = 1/1.1658e-11;

    GC_equation(Particle &p){
        q = p.q; m = p.m;
    }
 /*
    y[0-2]  =   GC_position
    y[3]    =   u
    y[4]    =   gyrophase
    y[5-7]  =   B_effective
    y[8-10] =   E_eff_cross_B_hat
    y[11]   =   B_eff_dot_E_eff
    y[12]   =   B_eff_parallell
    y[13]   =   B_eff_amplitude
 */
    void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydx){

        for(int i = 0; i < 3; i++) dydx[i] = (y[3]/y[12])*y[i+5] + (c/y[12])*y[i+8];            // dydx[0-2] = GC_v[0-2]

        dydx[3] = (q / m) * (1 / y[12]) * y[10];                                                // dydx[3] = d/dt u

        dydx[4] = (q / (m*c)) * y[13];                                                          // dydx[4] = d/dt gyrophase

        for(int i = 5; i<y.size();i++) dydx[i] = 0;                                             // Rest are set in integrate_GC()
    }

};




#endif  