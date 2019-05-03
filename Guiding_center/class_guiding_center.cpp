#include "Guiding_center/class_guiding_center.h"

#include "Functions/functions.h"

int Guiding_Center::initialize_new_GC(Particle &particle, Bfield &bfield, Ran &rng, double t){
      particle.initialize_new_particle(rng);
  //testing:
  particle.v = { 0.9*c, 0, sqrt(1-0.9*0.9)*c };
 
  bfield.generate_bfield_at_point(t, particle.pos);
 
  gyrofrequency = 8.987 * (particle.q / particle.m) * GCT::vector_amplitude(bfield.B);  // 8.987 is the unit coefficient for µG, e and MeV/c^2


 // Calculate initial GC-velocity
  v_parallell = calculate_v_parallell(particle.v, bfield.B_hat);  


  GC_velocity = { v_parallell*bfield.B_hat[0], v_parallell*bfield.B_hat[1], v_parallell*bfield.B_hat[2] };


  u = GCT::vector_dot_product(GC_velocity, bfield.B_hat);


 // Calculate initial GC-position
  v_perp = calculate_v_perp(particle.v, v_parallell);
  R_Larmor = 1.0810076e-15 * (v_perp / c) * (particle.E / bfield.B_amp_current());
 
  v_perp_hat = { particle.v[0] - GC_velocity[0], particle.v[1] - GC_velocity[1], particle.v[2] - GC_velocity[2] };
 
  GCT::normalize_vector(v_perp_hat);                                                  // This is now v_perp_hat

  a_hat = GCT::vector_cross_product(v_perp_hat, bfield.B_hat);
  GCT::normalize_vector(a_hat);

  GC_position = { particle.pos[0] + R_Larmor*a_hat[0], particle.pos[1] + R_Larmor*a_hat[1], particle.pos[2] + R_Larmor*a_hat[2] };


  hat_1 = { -a_hat[0], -a_hat[1], -a_hat[2] };

  return 0;
}

double Guiding_Center::calculate_v_parallell(const std::vector<double> &particle_velocity, const std::vector<double> &B_hat){
  return GCT::vector_dot_product(particle_velocity, B_hat);
}

double Guiding_Center::calculate_v_perp(const std::vector<double> &particle_velocity, const double v_parallell){
  return sqrt(pow(GCT::vector_amplitude(particle_velocity), 2) - pow(v_parallell, 2));
}



Guiding_Center::Guiding_Center(Initializer &init){
    GC_velocity.resize(3,0), GC_position.resize(3,0);     // Guiding center velocity and position
}