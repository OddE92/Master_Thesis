#include "Particle/class_particle.h"

#include <iostream>

int Particle::initialize_new_particle(Ran &rng){
    
    // Avoid using 0 due to singularity in some bfields at r = 0
    pos = {100, 100, 0};
    //pos = { 0, 0, 0, };

    double ranPhi = 2 * M_PI * rng.doub();
    double ranTheta = acos(1.0 - 2 * rng.doub());
    
    v[0] = -cos(ranPhi) * sin(ranTheta) * v_total;      // -signs becuse the 0-seeded RNG-velocity points inwards in the spiral field
    v[1] = -sin(ranPhi) * sin(ranTheta) * v_total;
    v[2] = -cos(ranTheta) * v_total;

    return 0;
}

Particle::Particle(Initializer &init){

    E = init.E;                             // [E] = [eV]
    q = init.q;                             // q = Z when Q = Z * e, [q] = [e]
    m = init.m;                             // [m] = [MeV/c^2]
    gamma_l = E/(1.0e6 * m);                // Lorentz factor calculated from E = gamma_l * m * c^2, with [m] = [MeV/c^2]

    //v_total = GCT::c * std::sqrt(1 - (1 / std::pow(gamma_l, 2)));     // Relativistic
    v_total = std::sqrt( (2*E*GCT::c*GCT::c) / (1.0e6*m));

    gamma_l_NR = 1/std::sqrt(1-std::pow(v_total/GCT::c, 2));
}

Particle::Particle(){

    // Default is a proton with E = 10^17eV

    std::cout << "Particle: Default constructor called. \n\n";
    
    E = 1.0e17;
    q = 1;                                  // Proton charge
    m = 938.2720813;                        // Proton mass
    gamma_l = E/(1.0e6 * m);                // Lorentz factor calculated from E = gamma_l * m * c^2, with [m] = [MeV/c^2]
    
    //Find total initial velocity
    v_total = GCT::c * std::sqrt(1 - (1 / std::pow(gamma_l, 2)));

}