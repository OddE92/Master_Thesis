#ifndef CLASS_PARTICLE
#define CLASS_PARTICLE

#define _USE_MATH_DEFINES

#include "Initializer/initializer.h"
#include "NR3/ran.h"

#include <vector>
#include <cmath>

constexpr double c = 2.99792458e8;              // Speed of light in m/s

class Particle{
    public:
        std::vector<double> pos, v;

        double E, q, m, gamma_l;
        double v_total;

        int initialize_new_particle(Ran &rng);

        Particle();
        Particle(Initializer &init);
        ~Particle(){};
    private:
};

#endif