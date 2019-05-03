/**********************************************

 * This program is a translation of the structure.f90 code written by Kristian Joten Andersen
 * for his thesis "Charged Particle Trajectories in the Local Superbubble".

 * The program was translated to C++ by Odd-Einar C. Nervik for his project thesis, under the guidance of
 * Michael Kachlerieß.
 
 * Some of the comments in the code are taken directly from the structure.f90-program.
 * The program follows the method shown in:

 *      Charged-particle motion in multidimensional magnetic-field turbulence
 *      Giacalone, J. and Jokipii, J. R.
 *      Department of Planetary Sciences, university of Arizona, Tucson AZ
 *      May 1994

 * The final calculation of the turbulent field follows the equation in:
 
 *      The transport of cosmic rays across a turbulent magnetic field
 *      Giacalone, J. and Jokipii, J. R.
 *      Department of Planetary Sciences, university of Arizona, Tucson AZ
 *      February 1999

 * Which has a different rotation and translation direction than the 1994 article.

 * The output of the program are two files. The first three numbers are the spatial coordinates x, y and z.
 * The last three numbers are Bx, By and Bz respectively.
 * *** The part where a box volume with a magnetic field has been deprecated. Now the field is generated at the provided point instead.
 
 * The second file is scaled for improved readability when plotted.

**********************************************/

#ifndef CLASS_BFIELD
#define CLASS_BFIELD

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#include "Particle/class_particle.h"

#include "NR3/ran.h"
#include "Initializer/initializer.h"

struct Bfield_func{
    
 /***** Insert the B-field function component-wise in x, y and z *****/

     const double theta = 11.5* M_PI / 180.0;
     double phi;

     double x(double t, std::vector<double> &pos){
        phi = std::atan2(pos[1], pos[0]);
        return (std::sin(theta) * std::cos(phi) - std::cos(theta) * std::sin(phi));
     }
     double y(double t, std::vector<double> &pos){
        phi = std::atan2(pos[1], pos[0]);
        return (std::sin(theta) * std::sin(phi) + std::cos(theta) * std::cos(phi));
     }
     double z(double t, std::vector<double> &pos){
         return 0;
     }

};


/********************************/
/********* CLASS BFIELD *********/

class Bfield{
   public:
      //std::array<double, 3> turbAtPoint, B, B_effective, E_effective;
      std::vector<double> turbAtPoint, B, B_hat;
      std::vector<double> B_effective, E_effective;
      std::vector<std::vector<double> > B_hat_partial;

      bool gen_turb = true;

      //Functions for generating the B-field
      int generate_bfield_at_point(double t, std::vector<double> &B_out, std::vector<double> &pos);
      int generate_bfield_at_point(double t, std::vector<double> &pos);

      int calculate_B_effective(std::vector<double> &GC_velocity, std::vector<double> &pos, double u);

      int calculate_E_effective(Particle &particle, double v_perp);

      int calculate_partial_b_hat(      std::vector<double> &B_at_point, std::vector<double> GC_velocity,
                                        std::vector<double> &GC_position, double timestep, double t);

      double B_amp_current();

      //General functions
      int generate_turbulence_at_point(std::vector<double> &pos);
      int initialize_turbulence(Ran &ran);
      int reinitialize_turbulence(Ran &ran);

      //Constructors and destructor
      Bfield();
      Bfield(Initializer &init, Ran &ran);
      ~Bfield();

 /********** CONSTANTS **********/

      //pi = M_PI
      const double two_pi = 2 * M_PI;

      std::complex<double> im; //For working with complex numbers

      //Magnetic field
      int n_k = 50; //Nr. of modes used to generate the turbulence
      double B_0 = 1.0;
      double B_rms_turb = 1.0;            //Normalized magnetic field,
      const double gamma = 5.0 / 3.0;     //Power law for the fluctuation spectrum
      double lambda_min = 0.2;            //Smallest wavelength, in parsecs
      double lambda_max = 10.0;           //Largest wavelength, in parsecs
      double k_min = two_pi / lambda_max; //Smallest wavenumber
      double k_max = two_pi / lambda_min; //Largest wavenumber

 /********* END CONSTANTS ********/

   private:

      Bfield_func B_static;

      //functions to initialize the turbulence. Only run these once!
      void initialize_phases(Ran &ran);
      void initialize_phases_from_file();
      void initialize_normalization();

      //Bools to check if functions need to be run or not.
      bool turbulence_is_initialized = false;

 /********* RANDOM PHASES ********/

      std::vector<double> a, b, p, t, s;          //a = alpha, b = beta, p = phi, t = theta, s = sign
      std::vector<double> ca, sa, cp, sp, ct, st; //c = cos, s = sin; ca = cos(a(n_k));

 /******* END RANDOM PHASES ******/

 /***** NORMALIZED PARAMETERS ****/

      std::vector<double> B_k, k;
      double dB_min;

 /*** END NORMALIZED PARAMETERS **/
};

#endif