#ifndef CLASS_GUIDING_CENTER
#define CLASS_GUIDING_CENTER

#include "Initializer/initializer.h"
#include "Particle/class_particle.h"
#include "Bfield/class_bfield.h"

#include "Units/units.h"

class Guiding_Center{
    public:
        std::array<Parsec, 3>   GC_position;
        std::array<mps, 3>      GC_velocity;
        std::array<double, 3>   v_perp_hat, a_hat;

        double gyrofrequency;
        double gyrophase;
        mps u;
        double mu;
    
        Parsec R_Larmor;
        Seconds timestep;
        _ps dudt;

        // Functions
        int initialize_new_GC(Particle &particle, Bfield &bfield, Ran &rng, Seconds t);

        double calculate_v_parallell(const std::array<mps, 3> &particle_velocity, const std::array<double, 3> &B_hat);
        double calculate_v_perp(const std::array<mps, 3> &particle_velocity, mps v_parallell);


        Guiding_Center(){};
        Guiding_Center(Initializer &init);
        ~Guiding_Center(){};

    private:
};






#endif