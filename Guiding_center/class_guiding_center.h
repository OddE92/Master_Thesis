#ifndef CLASS_GUIDING_CENTER
#define CLASS_GUIDING_CENTER

#include "Initializer/initializer.h"
#include "Particle/class_particle.h"
#include "Bfield/class_bfield.h"

class Guiding_Center{
    public:
        std::vector<double> GC_velocity, GC_position;
        std::vector<double> a_hat, hat_1;
        std::vector<double> v_perp_hat;

        double gyrofrequency;
        double gyrophase;
        double u;
        double v_perp, v_parallell;
        double R_Larmor;

        // Functions
        int initialize_new_GC(Particle &particle, Bfield &bfield, Ran &rng, double t);

        double calculate_v_parallell(const std::vector<double> &particle_velocity, const std::vector<double> &B_hat);
        double calculate_v_perp(const std::vector<double> &particle_velocity, const double v_parallell);


        Guiding_Center();
        Guiding_Center(Initializer &init);
        ~Guiding_Center(){};

    private:
};






#endif