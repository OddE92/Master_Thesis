#ifndef CLASS_PARTICLE
#define CLASS_PARTICLE

#define _USE_MATH_DEFINES

#include "Initializer/initializer.h"
#include "NR3/ran.h"
#include "Constants/constants.h"

#include "Units/units.h"

#include <vector>
#include <cmath>
#include <array>

class Particle{
    public:
        std::array<Parsec, 3>   pos;
        std::array<mps, 3>      v;

        double E, q, m, gamma_l;
        mps v_total, v_perp, v_parallell;
        double p_para;

        double gamma_l_NR;

        int initialize_new_particle(Ran &rng);

        Particle();
        Particle(Initializer &init);
        ~Particle(){};
    private:
};

#endif