#include "Guiding_center/class_guiding_center.h"

#include "Functions/functions.h"
#include "Trajectory/class_trajectory.h"

int Guiding_Center::initialize_new_GC(Particle &particle, Bfield &bfield, Ran &rng, double t){
  particle.initialize_new_particle(rng);              // Need a new particle to initiate a new GC
 
  std::array<double, 3> r = {particle.pos[0]*GCT::mtopc, particle.pos[1]*GCT::mtopc, particle.pos[2]*GCT::mtopc};

  bfield.generate_bfield_at_point(t, r);   // Generate the field at the starting point

 
  gyrofrequency = 8.987 * (particle.q / particle.m) * GCT::vector_amplitude(bfield.B);  // 8.987 is the unit coefficient for µG, e and MeV/c^2


 // Calculate initial GC-velocity
  v_parallell = calculate_v_parallell(particle.v, bfield.B_hat);          // v_parallell = v * \hat{b}

      // GC.v = \vect{v_parallell} = v_parallell * \hat{b} at start point
  GC_velocity = { v_parallell*bfield.B_hat[0], v_parallell*bfield.B_hat[1], v_parallell*bfield.B_hat[2] };


  u = GCT::vector_dot_product(GC_velocity, bfield.B_hat);                 // In general u = \vect{\dot{X}} * \hat{b}


 // Calculate initial GC-position
  v_perp = calculate_v_perp(particle.v, v_parallell);
  R_Larmor = 1.0810076e-15 * (v_perp / GCT::c) * (particle.E / GCT::vector_amplitude(bfield.B));
 
    // \vect{v_perp} = \vect{v} - \vect{v_parallell (=GC.v)}
  v_perp_hat = { particle.v[0] - GC_velocity[0], particle.v[1] - GC_velocity[1], particle.v[2] - GC_velocity[2] };
 
  GCT::normalize_vector(v_perp_hat);                                      // This is now v_perp_hat

  a_hat = GCT::vector_cross_product(v_perp_hat, bfield.B_hat);            // \hat{a} is the vector from the particle to the GC
  GCT::normalize_vector(a_hat);


    // the GC position is a distance R_larmor away from the particle position, in the direction of \hat{a}, unit = meter
  GC_position = { particle.pos[0] + (R_Larmor*a_hat[0])/GCT::mtopc, 
                  particle.pos[1] + (R_Larmor*a_hat[1])/GCT::mtopc, 
                  particle.pos[2] + (R_Larmor*a_hat[2])/GCT::mtopc 
                };


  hat_1 = { -a_hat[0], -a_hat[1], -a_hat[2] };                            // \hat{1} is fixed in this orientation wrt B pr definition

  return 0;
}

double Guiding_Center::calculate_v_parallell(const std::array<double, 3> &particle_velocity, const std::array<double, 3> &B_hat){
  return GCT::vector_dot_product(particle_velocity, B_hat);
}

double Guiding_Center::calculate_v_perp(const std::array<double, 3> &particle_velocity, const double v_parallell){
  return sqrt(pow(GCT::vector_amplitude(particle_velocity), 2) - pow(v_parallell, 2));
}



Guiding_Center::Guiding_Center(Initializer &init){
    timestep = init.dt;       // make sure timestep is set.
}