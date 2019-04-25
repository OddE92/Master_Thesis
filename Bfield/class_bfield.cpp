#include "Bfield/class_bfield.h"

/********** GENERATE FIELD AT POINT **********/
int Bfield::generate_bfield_at_point(double t, std::vector<double> &B_out, std::vector<double> &pos){

/* 
    This functions takes the position and time, generates the turbulence (without reinitializing it) 
    and calculates the regular field at the position and time. B_0 is used to scale the regular field.
*/
    
    if(!turbulence_is_initialized){ 
        initialize_turbulence();
    }
    if(gen_turb) {
        generate_turbulence_at_point(pos);
    }

    B_out[0] = B_0 * B_static.x(t, pos) + turbAtPoint[0];
    B_out[1] = B_0 * B_static.y(t, pos) + turbAtPoint[1];
    B_out[2] = B_0 * B_static.z(t, pos) + turbAtPoint[2];

    return 0;
}
/********** END GENERATE FIELD AT POINT **********/

int Bfield::generate_turbulence_at_point(std::vector<double> &pos){
    std::vector<double> e_k(n_k);
    
    std::vector< std::complex<double> > delta_B(3);
    std::vector< std::complex<double> > F(n_k), epsilon_x(n_k), epsilon_y(n_k), epsilon_z(n_k), B_x_k(n_k), B_y_k(n_k), B_z_k(n_k);

    // This function takes a position r, and generate a "random" turbulence to the B-field

    // Follow the articles calculation of delta_B (or delta_Omega, where delta_B = delta_Omega * mc/q)

    for(int i = 0; i < n_k; i++){
        using namespace std::complex_literals;

        // e_k is the exponent of delta_B
        e_k[i] = k[i] * (st[i]*cp[i]*pos[0] + st[i]*sp[i]*pos[1] + ct[i]*pos[2]) + b[i];
        //std::cout << e_k[i] << std::endl;

        // F = B(k)*e^i(kz+beta), B_k = B(k)
        F[i].real(B_k[i] * cos(e_k[i]));  F[i].imag(B_k[i] * sin(e_k[i]));

        // if(x == 0){
        //     std::cout << F[i] << endl;
        // }


        epsilon_x[i] = ca[i] * ct[i] * cp[i] - 1i * s[i] * sa[i] * sp[i];
        epsilon_y[i] = ca[i] * ct[i] * sp[i] + 1i * s[i] * sa[i] * cp[i];
        epsilon_z[i] = -ca[i] * st[i];

        // if(x == 0){
        //     std::cout << epsilon_x[i] << endl;
        // }

        //Calculate B(k) in the x-, y- and z-directions
        B_x_k[i] = F[i] * epsilon_x[i];                      // Complex number multiplication!
        B_y_k[i] = F[i] * epsilon_y[i];
        B_z_k[i] = F[i] * epsilon_z[i];

        // if(x == 0){
        //     std::cout << B_x_k[i] << endl; //<< ' ' << B_y_k[i] << ' ' << B_z_k[i] << endl;
        // }
    }

    for(int i = 0; i < delta_B.size(); i++){                 // Set delta_B to 0 before we start summing
        delta_B[i].imag(0.0); delta_B[i].real(0.0);
    }

        //Calculates the final delta_B in the x-, y- and z-directions.
    for(int i = 0; i < n_k; i++){
        delta_B[0] += B_x_k[i];
        delta_B[1] += B_y_k[i];
        delta_B[2] += B_z_k[i];
    }

    // if(x == 0){
    //     std::cout << delta_B[0] << ' ' << delta_B[1] << ' ' << delta_B[2] << endl;
    // }

// As we are concerned about the real value of the B-field we need to calculate the real part of delta_B

    for(int i = 0; i < delta_B.size(); i++){
        turbAtPoint[i] = real(delta_B[i]);
    }

// B = m*c*Omega/q, Omega = Omega_0 + dOmega. Check units!
    return 0;
}

/**** END GENERATE TURBULENCE ***/


/********** INITIALIZE **********/

int Bfield::initialize_turbulence(){
    if (!turbulence_is_initialized){
        initialize_phases();
        //initialize_phases_from_file();                  //For generating same result when testing
        initialize_normalization();

        turbulence_is_initialized = true;
        return 0;
    }
    std::cout << "Already initialized \n\n";
    return 1;
}

int Bfield::reinitialize_turbulence(){

    initialize_phases();
    initialize_normalization();

    turbulence_is_initialized = true;
    
    return 0;
}


/********************************/
/****** INITIALIZE PHASES *******/

void Bfield::initialize_phases(){
    double r;

    for(int i = 0; i < n_k; i++){                           // From the article:
        a[i] = two_pi * ran.doub();                         // a = alpha
        b[i] = two_pi * ran.doub();                         // b = beta
        p[i] = two_pi * ran.doub();                         // p = phi
        //t[i] = M_PI * ran.doub();
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
        a[i] = two_pi * ran_nums[j+0];                         // a = alpha
        b[i] = two_pi * ran_nums[j+1];                         // b = beta
        p[i] = two_pi * ran_nums[j+2];                         // p = phi
        t[i] = std::acos(1 - 2*ran_nums[j+3]);                 // t = theta with p(t) = sin(t)/2, 0<t<pi

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
    dlog_k = logk_diff/(n_k - 1);

    
    //The for loop spreads the values of k equally along the logarithmic scale.
    for(int i = 0; i < n_k; i++){
        k[i] = pow(10, log10(k_min) + i * dlog_k);
        //std::cout << k[i] << std::endl;
    }

  double sum = 0;
  for(int i = 0; i < n_k; i++){
      sum += pow(k[i]/k_min, -gamma);
  }
  
  dB_min = B_rms_turb * sqrt(2 / sum);
  for(int i = 0; i < n_k; i++){
      B_k[i] = dB_min * pow(k[i] / k_min, -gamma/2);
      //std::cout << B_k[i] << endl;
  }

  //Normalization is now done after the value given for B_rms_turb (default 1.0).  
}
/****** END INITIALIZE NORM *****/
/********************************/


Bfield::Bfield() : ran(15321){
    
    // See class declaration for std::vector descriptions.

    std::cout << "Default constructor called. \n\n";

    turbAtPoint.resize(3);

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);
}



Bfield::Bfield( Initializer &init) : B_0(init.B_0), B_rms_turb(init.B_rms_turb), ran(15321 + init.seed){
    
    // See class declaration for std::vector descriptions.
    n_k = init.n_k;
  
    turbAtPoint.resize(3);
    B.resize(3);

    gen_turb = init.generate_turbulence;

    lambda_max = init.lambda_max;
    lambda_min = init.lambda_min;

    k_min = two_pi / lambda_max;                     //Smallest wavenumber
    k_max = two_pi / lambda_min;                     //Largest wavenumber

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);
    initialize_turbulence();

}

Bfield::~Bfield(){
    //do nothing
}
