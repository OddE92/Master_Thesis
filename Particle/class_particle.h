#ifndef CLASS_PARTICLE
#define CLASS_PARTICLE

#define _USE_MATH_DEFINES

#include "Initializer/initializer.h"
#include "NR3/ran.h"
#include "Constants/constants.h"

#include <vector>
#include <cmath>

class Particle{
    public:
        std::array<double, 3> pos, v;

        double E, q, m, gamma_l;
        double v_total;

        int initialize_new_particle(Ran &rng);

        Particle();
        Particle(Initializer &init);
        ~Particle(){};
    private:
};

#endif