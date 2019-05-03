#ifndef INITIALIZER
#define INITIALIZER

#include <vector>

/******* INITIALIZING LIST ******/
struct Initializer{
    bool GCT;

   // Trajectory
    double t_start, t_end, dt;
    int t_end_y;
    int N;
    int num_points_recorded;

   // Magnetic field
    int n_k;
    double B_0, B_rms_turb;
    double lambda_max, lambda_min;
    double min_err, max_err;
    bool generate_turbulence;

   // Particle
    std::vector<double> start_pos;
    std::vector<double> start_vel;
    double E, q, m;
    
   // RNG
    int seed;
    int procID;
    
};
/***** END INITIALIZING LIST ****/


#endif