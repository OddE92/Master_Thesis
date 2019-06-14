#include "Bfield/class_bfield.h"

#include "Functions/functions.h"

/********** GENERATE FIELD AT POINT **********/

int Bfield::generate_bfield_at_point(Seconds t, std::array<microGauss, 3> &B_out, std::array<Parsec, 3> &pos){

/* 
    This functions takes the position and time, generates the turbulence (without reinitializing it) 
    and calculates the regular field at the position and time. B_0 can be used to scale the regular field.

    This version takes an external array to store the B-field.
*/
    
    if(!turbulence_is_initialized && gen_turb){ 
        std::cout << "Turbulence not initialized";
        throw("Turbulence not initialized");
    }
    if(gen_turb) {
        generate_turbulence_at_point(pos);
    }

    B_out[0] = B_0 * B_static.x(t, pos) + turbAtPoint[0];
    B_out[1] = B_0 * B_static.y(t, pos) + turbAtPoint[1];
    B_out[2] = B_0 * B_static.z(t, pos) + turbAtPoint[2];

    microGauss B_amp = GCT::vector_amplitude(B_out);

    if(B_amp == 0){
        std::cout << "In generate_bfield_at_point: B_amp == 0.\n";
        throw("B_amp = 0");
    }

    B_hat = { B[0]/B_amp, B[1]/B_amp, B[2]/B_amp };

    return 0;
}

// overload to set bfield.B instead of a passed vector
int Bfield::generate_bfield_at_point(Seconds t, std::array<Parsec, 3> &pos){

/* 
    This functions takes the position and time, generates the turbulence (without reinitializing it) 
    and calculates the regular field at the position and time. B_0 can be used to scale the regular field.

    This function uses Bfield.B to store the generated magnetic field.
*/
    
    if(!turbulence_is_initialized && gen_turb){ 
        std::cout << "Turbulence not initialized";
        throw("Turbulence not initialized");
    }
    if(gen_turb) {
        generate_turbulence_at_point(pos);
    }

    B[0] = B_0 * B_static.x(t, pos) + turbAtPoint[0];
    B[1] = B_0 * B_static.y(t, pos) + turbAtPoint[1];
    B[2] = B_0 * B_static.z(t, pos) + turbAtPoint[2];

    microGauss B_amp = GCT::vector_amplitude(B);

    if(B_amp == 0){
        std::cout << "In generate_bfield_at_point: B_amp == 0.\n";
        throw("B_amp = 0");
    }


    B_hat = { B[0]/B_amp, B[1]/B_amp, B[2]/B_amp };

    return 0;
}
/********** END GENERATE FIELD AT POINT **********/

int Bfield::generate_turbulence_at_point(std::array<Parsec, 3> &pos){

    // This function takes a position r, and generate a "random" turbulence to the B-field

    // Follow Jokipii 1994s calculation of delta_B (or delta_Omega, where delta_B = delta_Omega * mc/q)

    for(int i = 0; i < n_k; i++){
        using namespace std::complex_literals;

        // e_k is the exponent of delta_B
        e_k[i] = k[i] * (st[i]*cp[i]*pos[0] + st[i]*sp[i]*pos[1] + ct[i]*pos[2]) + b[i];


        // F = B(k)*e^i(kz+beta), B_k = B(k)
        F[i].real(B_k[i] * cos(e_k[i]));  F[i].imag(B_k[i] * sin(e_k[i]));

        // epsilon = cos(alpha)\hat{x'} + i*sin(alpha)\hat{y'}
        epsilon_x[i] = ca[i] * ct[i] * cp[i] - 1i * s[i] * sa[i] * sp[i];
        epsilon_y[i] = ca[i] * ct[i] * sp[i] + 1i * s[i] * sa[i] * cp[i];
        epsilon_z[i] = -ca[i] * st[i];


        //Calculate B(k) in the x-, y- and z-directions
        B_x_k[i] = F[i] * epsilon_x[i];                      // Complex number multiplication!
        B_y_k[i] = F[i] * epsilon_y[i];
        B_z_k[i] = F[i] * epsilon_z[i];

    } //End for

    for(int i = 0; i < delta_B.size(); i++){                 // Set delta_B to 0 before we start summing
        delta_B[i].imag(0.0); delta_B[i].real(0.0);
    }

        //Calculates the final delta_B in the x-, y- and z-directions.
    for(int i = 0; i < n_k; i++){
        delta_B[0] += B_x_k[i];
        delta_B[1] += B_y_k[i];
        delta_B[2] += B_z_k[i];
    }


        // As we are concerned about the real value of the B-field we need to calculate the real part of delta_B
 
    for(int i = 0; i < delta_B.size(); i++){
        turbAtPoint[i] = real(delta_B[i]);
    }

        // B = m*c*Omega/q. Check units!
    return 0;
}

/**** END GENERATE TURBULENCE ***/

/* CALCULATE PARTIAL DERIVATIVE */
int Bfield::calculate_partial_b_hat(    std::array<microGauss, 3> &B_at_point, std::array<mps, 3> &GC_velocity, 
                                        std::array<Parsec, 3> &GC_position, Seconds t, Seconds timestep){

    // B_hat_partial[0][i] = d/di B__hat_x, [1][i] = d/di B__hat_y, [2][i] = d/di B__hat_z

    std::array<Parsec, 3> GC_next_position;                      // Holds GC_position after taking a step in x, y or z-direction
    std::array<microGauss, 3> B_at_next_point;

    Meter step_size = GCT::vector_amplitude(GC_velocity) * timestep; 
  

 // i loops through dx, dy and dz
    for(int i = 0; i < GC_next_position.size(); i++){                       // Step through x, y, z as [0], [1], [2].

        for (int j = 0; j < GC_next_position.size(); j++){          

            if(i == j)
                GC_next_position[j] = GC_position[j] + step_size*GCT::mtopc;           // Take a step in the i-direction only
            else
                GC_next_position[j] = GC_position[j];                   
        
        } //End for j
     
    
        generate_bfield_at_point(t, B_at_next_point, GC_next_position);     // Generate field at new position
     
      // This is needed in calc_E*, but we save it here to avoid having to generate the B-field 3 extra times.
        Nabla_B[i] = units::mGtoT*(GCT::vector_amplitude(B_at_next_point) - GCT::vector_amplitude(B_at_point)) / step_size;

        microGauss B_norm_p  = GCT::vector_amplitude(B_at_point);
        microGauss B_norm_np = GCT::vector_amplitude(B_at_next_point);
     

     // Calculate d/di(b_hat) = b_hat(r+dr) - b_hat(r) / (2*dr), with dr = dx, dy, dz for [0], [1], [2]. 
        for(int j = 0; j < B_at_next_point.size(); j++){
            B_hat_partial[i][j] = (B_at_next_point[j]/B_norm_np - B_at_point[j]/B_norm_p) / step_size;
            
        } //End for j

    }//End for i
    
    return 0;
}

/* END CALC PARTIAL DERIVATIVE  */

/** CALCULATE EFFECTIVE FIELDS **/

int Bfield::calculate_B_effective(std::array<mps, 3> &GC_velocity, std::array<Parsec, 3> &pos, mps u, MeVc2 m, e_charge q){

    // Calculates B* using the identity 
    // (Nabla cross scalar*vect) = scalar*(Nabla cross vect) + (Nabla*u) cross vect

    // B* = B + (m*c)/q * (Nabla cross u*\hat{b}) 
 //std::cout << "\nIn calc_B*:\nGC_velocity: " << GC_velocity << std::endl;
 //std::cout << "B_hat_partial:\n" << B_hat_partial[0] << '\n' << B_hat_partial[1] << '\n' << B_hat_partial[2] << std::endl;

    
 
        // Nabla*u = \dot{X} * nabla * \hat{b}
    for(int i = 0; i < nabla_u.size(); i++){
        nabla_u[i] = GCT::vector_dot_product(GC_velocity, B_hat_partial[i]);
    }
 //std::cout << "Nabla_u: " << nabla_u << std::endl;

    std::array<_ps, 3> nabla_u_cross_b_hat = GCT::vector_cross_product(nabla_u, B_hat);
 //std::cout << "Nabla_u_cross_b_hat: " << nabla_u_cross_b_hat << std::endl;
    std::array<double, 3> nabla_cross_b_hat = {
        B_hat_partial[1][2] - B_hat_partial[2][1],
        B_hat_partial[2][0] - B_hat_partial[0][2],
        B_hat_partial[0][1] - B_hat_partial[1][0]
    };
 //std::cout << "Nabla_cross_b_hat: " << nabla_cross_b_hat << std::endl;
 //std::cout << "u*Nabla_cross_b_hat: " << u*nabla_cross_b_hat[0] << ' ' << u*nabla_cross_b_hat[1] << ' ' << u*nabla_cross_b_hat[2] << std::endl;
        // B* = B + (m*c)/q * (Nabla cross u*\hat{b})               GCT::mq and GCT::c are defined in constants.h
    B_effective = {
        B[0] + (1e10*GCT::mq)*(m/q)*(u*nabla_cross_b_hat[0] + nabla_u_cross_b_hat[0]),
        B[1] + (1e10*GCT::mq)*(m/q)*(u*nabla_cross_b_hat[0] + nabla_u_cross_b_hat[1]),
        B[2] + (1e10*GCT::mq)*(m/q)*(u*nabla_cross_b_hat[0] + nabla_u_cross_b_hat[2])
    };
 //std::cout << "Parenthesis in B*: \n " << (1e10*GCT::mq)*(m/q)*(u*nabla_cross_b_hat[0] + nabla_u_cross_b_hat[0]) << std::endl;
 //std::cout << ' ' << (1e10*GCT::mq)*(m/q)*(u*nabla_cross_b_hat[1] + nabla_u_cross_b_hat[1]) << std::endl;
 //std::cout << ' ' << (1e10*GCT::mq)*(m/q)*(u*nabla_cross_b_hat[2] + nabla_u_cross_b_hat[2]) << std::endl;

    return 0;
}

int Bfield::calculate_E_effective(  const Coulomb &q, const double mu){

    // Calculates E* using
    // E* = - (µ/q) * (Nabla * B),  with µ = (m * v_perp^2) / (2*B) := const

    // Possibly needs to be changed to 
    // E* = -(m/2q)* (Nabla * v_perp^2)

    // Nabla_B calculated in B_hat_partial for efficiency.


    E_effective = {
            - (mu / q)  * Nabla_B[0],
            - (mu / q)  * Nabla_B[1],
            - (mu / q)  * Nabla_B[2]
        };

    // testing <B> = const:
    //E_effective = { 0, 0, 0 };
    return 0;
}


int Bfield::calculate_E_effective(const Kg &m, const Coulomb &q, mps v_perp){
    

    // v_perp = v - v_parallel, v_parallel = v*\hat{b}      v is assumed const
    // => Nabla * v_perp = - v * Nabla \hat{b} = - v_part * nabla_u

    E_effective = {
        (m / q) * v_perp * nabla_u[0],
        (m / q) * v_perp * nabla_u[1],
        (m / q) * v_perp * nabla_u[2]
    };

    return 0;
}

/** END CALC EFFECTIVE FIELDS  **/


/********** INITIALIZE **********/

int Bfield::initialize_turbulence(Ran &ran){
    if (!turbulence_is_initialized){
        initialize_phases(ran);
        //initialize_phases_from_file();                    //For generating same result when testing
        initialize_normalization();

        turbulence_is_initialized = true;
        return 0;
    }
    std::cout << "Already initialized \n\n";
    return 1;
}

int Bfield::reinitialize_turbulence(Ran &ran){              // Forcing new initialization if used.

    initialize_phases(ran);
    initialize_normalization();

    turbulence_is_initialized = true;
    
    return 0;
}


/********************************/
/****** INITIALIZE PHASES *******/

void Bfield::initialize_phases(Ran &ran){
    double r;

    for(int i = 0; i < n_k; i++){                           // From the article:
        a[i] = two_pi * ran.doub();                         // a = alpha
        b[i] = two_pi * ran.doub();                         // b = beta
        p[i] = two_pi * ran.doub();                         // p = phi
        t[i] = std::acos(1 - 2*ran.doub());                 // t = theta with p(t) = sin(t)/2, 0<t<pi

        s[i] = 1.0;                                         // s = sign in epsilon (+-)
        r = ran.doub();
        if(r > 0.5) s[i] = -1.0;

        cp[i] = cos(p[i]);                                  //sines and cosines of theta, phi and alpha
        sp[i] = sin(p[i]);
        ct[i] = cos(t[i]);
        st[i] = sin(t[i]);
        ca[i] = cos(a[i]);
        sa[i] = sin(a[i]);
    }
}

void Bfield::initialize_phases_from_file(){
    using namespace std;
    int j = 0;
    double r;
    string str;
    std::vector<double> ran_nums(10000);

    ifstream file;

    file.open("RNG.txt");

    while(getline(file, str)){
        ran_nums[j] = stod(str);
        j++;
    }

    file.close();

    j=0;

    for(int i = 0; i < n_k; i++){                           // From the article:
        a[i] = two_pi * ran_nums[j+0];                      // a = alpha
        b[i] = two_pi * ran_nums[j+1];                      // b = beta
        p[i] = two_pi * ran_nums[j+2];                      // p = phi
        t[i] = std::acos(1 - 2*ran_nums[j+3]);              // t = theta with p(t) = sin(t)/2, 0<t<pi

        s[i] = 1.0;                                         // s = sign in epsilon (+-)
        r = ran_nums[j+4];
        if(r > 0.5) s[i] = -1.0;

        cp[i] = cos(p[i]);                                  //sines and cosines of theta, phi and alpha
        sp[i] = sin(p[i]);
        ct[i] = cos(t[i]);
        st[i] = sin(t[i]);
        ca[i] = cos(a[i]);
        sa[i] = sin(a[i]);

        j+=5;
    }

}
/***** END INITIALIZE PHASES ****/
/********************************/

/********************************/
/*** INITIALIZE NORMALIZATION ***/

void Bfield::initialize_normalization(){

    double logk_diff, dlog_k;

    logk_diff = log10(k_max) - log10(k_min);
    dlog_k = logk_diff/(n_k - 1);                       // Find the logarithmic step-size for k

    
    //The for loop spreads the values of k equally along the logarithmic scale.
    for(int i = 0; i < n_k; i++){
        k[i] = pow(10, log10(k_min) + i * dlog_k);
        //std::cout << k[i] << std::endl;
    }

  double sum = 0;                                       // calculate sum((k_i / k_min)^-gamma)
  for(int i = 0; i < n_k; i++){
      sum += pow(k[i]/k_min, -gamma);
  }
  
  dB_min = B_rms_turb * sqrt(2 / sum);                  // Calculate dB_min from given RMS-value and the sum

  for(int i = 0; i < n_k; i++){                         // Calculate all B_k values
      B_k[i] = dB_min * pow(k[i] / k_min, -gamma/2);
  }

  //Normalization is now done after the value given for B_rms_turb (default 1.0).  
}
/****** END INITIALIZE NORM *****/
/********************************/


Bfield::Bfield(){

    std::cout << "Default constructor called. \n\n";

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);

    e_k.resize(n_k);
    
    F.resize(n_k); 
    epsilon_x.resize(n_k); epsilon_y.resize(n_k); epsilon_z.resize(n_k);
    B_x_k.resize(n_k); B_y_k.resize(n_k); B_z_k.resize(n_k);
}



Bfield::Bfield( Initializer &init, Ran &ran) : B_0(init.B_0), B_rms_turb(init.B_rms_turb){
    
    n_k = init.n_k;
 
    gen_turb = init.generate_turbulence;

    lambda_max = init.lambda_max;
    lambda_min = init.lambda_min;

    k_min = two_pi / lambda_max;                    // Smallest wavenumber
    k_max = two_pi / lambda_min;                    // Largest wavenumber

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);

    e_k.resize(n_k);

    F.resize(n_k); 
    epsilon_x.resize(n_k); epsilon_y.resize(n_k); epsilon_z.resize(n_k);
    B_x_k.resize(n_k); B_y_k.resize(n_k); B_z_k.resize(n_k);

    if(gen_turb){
        initialize_turbulence(ran);
    } else{
        B_rms_turb = 0.0;                           // Set to 0 for proper R_larmor-calculations
    }
}

/* This is generated by the compiler?
// Copy-constructor
Bfield::Bfield(const Bfield &b){
    turbAtPoint = b.turbAtPoint;
    B = b.B;
    B_hat = b.B_hat;
    B_effective = b.B_effective;
    E_effective = b.E_effective;
    B_hat_partial = b.B_hat_partial;

    gen_turb = b.gen_turb;

    turbulence_is_initialized = b.turbulence_is_initialized;

    n_k = b.n_k; 
    B_0 = b.B_0;
    B_rms_turb = b.B_rms_turb;            
    lambda_min = b.lambda_min;          
    lambda_max = b.lambda_max;           
    k_min = b.k_min; 
    k_max = b.k_max; 

    a = b.a; this->b = b.b; p = b.p; t = b.t; s = b.s; 
    ca = b.ca; sa = b.sa; cp = b.cp; sp = b.sp; ct = b.ct; st = b.st; 

    B_k = b.B_k; k = b.k;
    dB_min = b.dB_min;
}*/

Bfield::~Bfield(){
    //do nothing
}
