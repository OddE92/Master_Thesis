#ifndef ODEINT
#define ODEINT

#include "NR3/nr3.h"
#include "Bfield/class_bfield.h"
#include "Trajectory/class_trajectory.h"
#include "Initializer/initializer.h"
#include "Functions/functions.h"
#include "Guiding_center/class_guiding_center.h"

#include <iomanip>
#include <iostream>

/* 
    Only changes here are in the integrate functions, to store D_ij at each point.
*/

struct Output
{
    Int kmax;
    Int nvar;
    Int nsave;
    bool dense;
    Int count;
    Doub x1, x2, xout, dxout;
    VecDoub xsave;
    MatDoub ysave;
    Output() : kmax(-1), dense(false), count(0) {}
    Output(const Int nsavee) : kmax(500), nsave(nsavee), count(0), xsave(kmax)
    {
        dense = nsave > 0 ? true : false;
    }
    void init(const Int neqn, const Doub xlo, const Doub xhi)
    {
        nvar = neqn;
        if (kmax == -1)
            return;
        ysave.resize(nvar, kmax);
        if (dense)
        {
            x1 = xlo;
            x2 = xhi;
            xout = x1;
            dxout = (x2 - x1) / nsave;
        }
    }
    void resize()
    {
        Int kold = kmax;
        kmax *= 2;
        VecDoub tempvec(xsave);
        xsave.resize(kmax);
        for (Int k = 0; k < kold; k++)
            xsave[k] = tempvec[k];
        MatDoub tempmat(ysave);
        ysave.resize(nvar, kmax);
        for (Int i = 0; i < nvar; i++)
            for (Int k = 0; k < kold; k++)
                ysave[i][k] = tempmat[i][k];
    }
    template <class Stepper>
    void save_dense(Stepper &s, const Doub xout, const Doub h)
    {
        if (count == kmax)
            resize();
        for (Int i = 0; i < nvar; i++)
            ysave[i][count] = s.dense_out(i, xout, h);
        xsave[count++] = xout;
    }
    void save(const Doub x, VecDoub_I &y)
    {
        if (kmax <= 0)
            return;
        if (count == kmax)
            resize();
        for (Int i = 0; i < nvar; i++)
            ysave[i][count] = y[i];
        xsave[count++] = x;
    }
    template <class Stepper>
    void out(const Int nstp, const Doub x, VecDoub_I &y, Stepper &s, const Doub h)
    {
        if (!dense){
            std::cout << "Dense output not set in Output!" << '\n';
            throw("dense output not set in Output!");
        }
        if (nstp == -1)
        {
            save(x, y);
            xout += dxout;
        }
        else
        {
            while ((x - xout) * (x2 - x1) > 0.0)
            {
                save_dense(s, xout, h);
                xout += dxout;
            }
        }
    }
};




template <class Stepper>
struct Odeint
{
    static const Int MAXSTP = 500000;

    Doub EPS;
    Int nok;
    Int nbad;
    Int nvar;
    Doub x1, x2, hmin;
    bool dense;
    VecDoub y, dydx;
    VecDoub &ystart;
    Output &out;
    typename Stepper::Dtype &derivs;
    Stepper s;
    Int nstp;
    Doub x, h;
    Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
           const Doub atol, const Doub rtol, const Doub h1,
           const Doub hminn, Output &outt, typename Stepper::Dtype &derivss, 
           Bfield &bfield, Particle &particle
           );
    Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
           const Doub atol, const Doub rtol, const Doub h1,
           const Doub hminn, Output &outt, typename Stepper::Dtype &derivss);
    void integrate();
    void integrate( std::vector<double> &D_ij, std::vector<double> &D_ij_time, std::vector<double> &r_vect, 
                    Particle &particle, const double x_max);
    void integrate( std::vector<double> &D_ij, std::vector<double> &D_ij_time, std::vector<double> &r_vect, 
                    Particle &particle, Bfield &bfield, const double x_max);                
    void integrate_GC(  Trajectory &trajectory, const double x_max);            // GC passed through GC_equation instead
};

// Original constructor below:
template <class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
                        const Doub atol, const Doub rtol, const Doub h1, const Doub hminn,
                        Output &outt, typename Stepper::Dtype &derivss) :  nvar(ystartt.size()),
                                                                          y(nvar), dydx(nvar), ystart(ystartt), x(xx1), nok(0), nbad(0),
                                                                          x1(xx1), x2(xx2), hmin(hminn), dense(outt.dense), out(outt), derivs(derivss),
                                                                          s(y, dydx, x, atol, rtol, dense)
{
    EPS = std::numeric_limits<Doub>::epsilon();
    h = SIGN(h1, x2 - x1);
    for (Int i = 0; i < nvar; i++)
        y[i] = ystart[i];
    out.init(s.neqn, x1, x2);
}

/************************************************************/
/******************** INTEGRATE FUNCTION ********************/

// After the last changes integrate_GC and integrate are almost identical.

template <class Stepper>
void Odeint<Stepper>::integrate_GC(Trajectory &trajectory, const double x_max)         //x = t for this program
{
    int count = 0, logcount = 1;                            // Counters to track when to save D_ij
    int posCounter = 0;

    double t;                                               // Used to test if step should be recorded
    //std::ofstream file;
    //file.open("Data/GC_trajectory_ODE.dat");
    
 //std::cout << "y[0-2] at start: " << y[0] << ' ' << y[1] << ' ' << y[2] << std::endl;

    //derivs.r = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
    //derivs.GC.GC_velocity = { y[3]*derivs.bfield.B_hat[0], y[3]*derivs.bfield.B_hat[1], y[3]*derivs.bfield.B_hat[2] };
      
    //derivs.bfield.calculate_partial_b_hat(derivs.bfield.B, derivs.GC.GC_velocity, derivs.r, x, derivs.GC.timestep);

    derivs(x, y, dydx);                                     // First step

    derivs.GC.timestep = h; 
    
    if (dense)
        out.out(-1, x, y, s, h);
    else
        out.save(x, y);

    
    for (nstp = 0; x <= x_max; nstp++)
    {

        if ((x + h * 1.0001 - x2) * (x2 - x1) > 0.0){
            h = x2 - x;
        }
     //std::cout << "\n\n/**********************************************************************************************/\n";

      //derivs.r = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
      //derivs.GC.GC_velocity = { y[3]*derivs.bfield.B_hat[0], y[3]*derivs.bfield.B_hat[1], y[3]*derivs.bfield.B_hat[2] };
      
      //derivs.bfield.calculate_partial_b_hat(derivs.bfield.B, derivs.GC.GC_velocity, derivs.r, x, derivs.GC.timestep);

      //derivs.bfield.calculate_E_effective(derivs.particle, derivs.GC.GC_velocity, derivs.GC.GC_position, derivs.GC.mu, derivs.GC.timestep, derivs.GC.dudt); 
      
      
      // Take a step
        s.step(h, derivs);


        if (s.hdid == h)
            ++nok;
        else{
            ++nbad;
        }
        if (dense)
            out.out(nstp, x, y, s, s.hdid);
        else{
          /*  std::cout << "******* SAVING Y ********\n";
            std::cout << "y: " << y[0] << ' ' << y[1] << ' ' << y[2] << ' ' << y[3] << ' ' << y[4] << std::endl;
            file << y[0]*GCT::mtopc << ' ' << y[1]*GCT::mtopc << ' ' << y[2]*GCT::mtopc << ' ' << y[3] << ' ' << y[4] << std::endl;
            std::cout << "dydx: " << s.dydx[0] << ' ' << s.dydx[1] << ' ' << s.dydx[2] << ' ' << s.dydx[3] << ' ' 
                      << s.dydx[4] << std::endl;
            std::cout << "GC.u: " << derivs.GC.u << std::endl;
            std::cout << "\nB_hat_partial:\n " << derivs.bfield.B_hat_partial[0] << "\n " <<  derivs.bfield.B_hat_partial[1] << "\n " << derivs.bfield.B_hat_partial[2] << std::endl;
            std::cout << "B: " << derivs.bfield.B << std::endl;
            std::cout << "B*: " << derivs.bfield.B_effective << std::endl;
            std::cout << "E*: " << derivs.bfield.E_effective << std::endl;
          */  out.save(x, y);
        }

        derivs.GC.timestep = h;                             // Make sure timestep is updated (GC_equation r.GC.timestep)

        derivs.GC.R_Larmor = 1.0810076e-15 * (derivs.GC.v_perp / GCT::c) * (derivs.particle.E / GCT::vector_amplitude(derivs.bfield.B));
        derivs.GC.GC_position = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
        //std::cout << "derivs timestep: " << derivs.GC.timestep << std::endl;
        

        if(std::isnan(y[0])) return;                        // Break if the position goes to NaN
        if(h < 1e-04) {std::cout << "h < 1e-4\n"; return;}

            // If t > count * 10^yearcount, store value of D_ij
        double t = x+h;
        if( t/31557600 - (logcount * pow(10, count)) > __DBL_EPSILON__){
            using GCT::mtopc;
            // Add each particles coordinate to the sum in D_ij
            // D_ij holds the upper half, including the diagonal, of a symmetric matrix
            // stored as a 1D-array to comply with MPIs send and recieve functions.
            trajectory.D_ij[posCounter + 0] += (y[0]*mtopc) * (y[0]*mtopc);                             //D_11
            trajectory.D_ij[posCounter + 1] += (y[0]*mtopc) * (y[1]*mtopc);                             //D_12
            trajectory.D_ij[posCounter + 2] += (y[0]*mtopc) * (y[2]*mtopc);                             //D_13
            trajectory.D_ij[posCounter + 3] += (y[1]*mtopc) * (y[1]*mtopc);                             //D_22
            trajectory.D_ij[posCounter + 4] += (y[1]*mtopc) * (y[2]*mtopc);                             //D_23
            trajectory.D_ij[posCounter + 5] += (y[2]*mtopc) * (y[2]*mtopc);                             //D_33
            trajectory.D_ij[posCounter + 6] += (y[0]*mtopc) * (y[0]*mtopc) + (y[1]*mtopc) * (y[1]*mtopc) + \
                                    (y[2]*mtopc) * (y[2]*mtopc); //D_average

            trajectory.D_ij_time[posCounter / 7] = x/31557600;

            trajectory.r_vect[3 * (posCounter / 7) + 0] = y[0]*mtopc;
            trajectory.r_vect[3 * (posCounter / 7) + 1] = y[1]*mtopc;
            trajectory.r_vect[3 * (posCounter / 7) + 2] = y[2]*mtopc;

            if(logcount == 9){
                count++;
                logcount = 1;
            }else{
                logcount++;
            }

            posCounter += 7;
        }//end if save D_ij




        if ((x - x2) * (x2 - x1) >= 0.0) //check if end is reached
        {   
            //std::cout << "I'm in this if thingy" << std::endl;
            //file.close();
            for (Int i = 0; i < nvar; i++)
                ystart[i] = y[i];
            if (out.kmax > 0 && abs(out.xsave[out.count - 1] - x2) > 100.0 * abs(x2) * EPS)
                out.save(x, y);

            return;
        }


     //std::cout << "\n/**********************************************************************************************/\n";

        if (abs(s.hnext) <= hmin){
            std::cout << "Step size too small in Odeint! \n";
            return;
        }
        h = s.hnext;
    }
   
    /*
    std::cout << "Too many steps! \n";
    throw("Too many steps in routine Odeint");
    */
}

template <class Stepper>
void Odeint<Stepper>::integrate(    std::vector<double> &D_ij, std::vector<double> &D_ij_time, std::vector<double> &r_vect,
                                    Particle &particle, const double x_max)                     //x = t for this program
{
    int count = 0, logcount = 1;                            // Counters to track when to save D_ij
    int posCounter = 0;

    double t;                                               // Used to test if step should be recorded

        // For testing E* and B*
    //std::array<double, 3> B, r;
    // Guiding_Center GC;

    //std::ofstream file;
    //file.open("Data/Test_E_B_eff.dat");   
    
    derivs(x, y, dydx);
    

    if (dense)
        out.out(-1, x, y, s, h);
    else
        out.save(x, y);

 std::cout << "odeint, x: " << x << std::endl;
 std::cout << "odeint, x_max: " << x_max << std::endl;
    for (nstp = 0; x <= x_max; nstp++)
    {

        if ((x + h * 1.0001 - x2) * (x2 - x1) > 0.0){
            h = x2 - x;
        }

      // Take a step
        s.step(h, derivs);


        if (s.hdid == h){
            ++nok;
        } else{
            ++nbad;
        }
        if (dense){
            out.out(nstp, x, y, s, s.hdid);
        } else{
            out.save(x, y);
        }


        /********** TEST B* AND E* ***********/
     /*   
        particle.pos = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
        particle.v = { y[3], y[4], y[5] };
     //std::cout << "\nparticle.v: " << particle.v << std::endl;
     //std::cout << "|.v|: " << GCT::vector_amplitude(particle.v) << std::endl;
        derivs.bfield.generate_bfield_at_point(0.0, particle.pos);
    
     //std::cout << "derivs.bfield.B: " << derivs.bfield.B << std::endl;
     //std::cout << "derivs.bfield.B_hat: " << derivs.bfield.B_hat << std::endl;
     //std::cout << "v * B_hat: " << y[3] - y[3] * derivs.bfield.B_hat[0] << ' ' <<
     //                           y[4] - y[4] * derivs.bfield.B_hat[1] << ' ' << y[5] - y[5] * derivs.bfield.B_hat[2] << std::endl;
        double dotprod = GCT::vector_dot_product(particle.v, derivs.bfield.B_hat);      
        std::array<double, 3> v_hat_perp =  {   y[3] - dotprod * derivs.bfield.B_hat[0],
                                                y[4] - dotprod * derivs.bfield.B_hat[1],
                                                y[5] - dotprod * derivs.bfield.B_hat[2]
                                            };

     //std::cout << "v_perp: " << v_hat_perp << std::endl;
     //std::cout << "|v_perp|: " << GCT::vector_amplitude(v_hat_perp) << std::endl;        
        std::array<double, 3> a_hat = GCT::vector_cross_product(v_hat_perp, derivs.bfield.B_hat);
        GCT::normalize_vector(a_hat);
     //std::cout << "a_hat: " << a_hat << std::endl;
        GC.R_Larmor = GCT::calculate_R_larmor(GCT::vector_amplitude(v_hat_perp), particle.E, GCT::vector_amplitude(derivs.bfield.B) );
     std::cout << "R_larmor: " << GC.R_Larmor << std::endl;

        GC.GC_position = {
            (y[0]*GCT::mtopc + a_hat[0]*GC.R_Larmor),
            (y[1]*GCT::mtopc + a_hat[1]*GC.R_Larmor),
            (y[2]*GCT::mtopc + a_hat[2]*GC.R_Larmor)
        };
        
        double v_para = GCT::vector_dot_product(particle.v, derivs.bfield.B_hat);
        GC.GC_velocity = {
            v_para * derivs.bfield.B_hat[0],
            v_para * derivs.bfield.B_hat[1],
            v_para * derivs.bfield.B_hat[2]
        };
     //std::cout << "v_para: " << v_para << std::endl;

     //std::cout << "sqrt v_para^2 + v_perp^2: " << std::sqrt(std::pow(v_para, 2) + std::pow(GCT::vector_amplitude(v_hat_perp), 2)) << std::endl;
        

        derivs.bfield.calculate_partial_b_hat(      derivs.bfield.B, GC.GC_velocity,
                                        GC.GC_position, 0.0, h);

        derivs.bfield.calculate_E_effective(particle, GCT::vector_amplitude(v_hat_perp));
        derivs.bfield.calculate_B_effective(GC.GC_velocity, GC.GC_position, GCT::vector_amplitude(GC.GC_velocity), particle.m, particle.q);
     //std::cout << "B*: " << derivs.bfield.B_effective << std::endl;
     //std::cout << "E*: " << derivs.bfield.E_effective << std::endl;

        file << GC.GC_position << ' ' << derivs.bfield.B_effective << ' ' << derivs.bfield.E_effective << std::endl;

     */
     


     //std::cout << "Timestep h: " << h << std::endl;
        double t = x+h;
        if( t/31557600 - (logcount * pow(10, count)) > __DBL_EPSILON__){
            using GCT::mtopc;
            // Add each particles coordinate to the sum in D_ij
            // D_ij holds the upper half, including the diagonal, of a symmetric matrix
            // stored as a 1D-array to comply with MPIs send and recieve functions.
            D_ij[posCounter + 0] += (y[0]*mtopc) * (y[0]*mtopc);                             //D_11
            D_ij[posCounter + 1] += (y[0]*mtopc) * (y[1]*mtopc);                             //D_12
            D_ij[posCounter + 2] += (y[0]*mtopc) * (y[2]*mtopc);                             //D_13
            D_ij[posCounter + 3] += (y[1]*mtopc) * (y[1]*mtopc);                             //D_22
            D_ij[posCounter + 4] += (y[1]*mtopc) * (y[2]*mtopc);                             //D_23
            D_ij[posCounter + 5] += (y[2]*mtopc) * (y[2]*mtopc);                             //D_33
            D_ij[posCounter + 6] += (y[0]*mtopc) * (y[0]*mtopc) + (y[1]*mtopc) * (y[1]*mtopc) + 
                                    (y[2]*mtopc) * (y[2]*mtopc);                             //D_average

            D_ij_time[posCounter / 7] = x/31557600;

            r_vect[3 * (posCounter / 7) + 0] = y[0]*mtopc;
            r_vect[3 * (posCounter / 7) + 1] = y[1]*mtopc;
            r_vect[3 * (posCounter / 7) + 2] = y[2]*mtopc;

            if(logcount == 9){
                count++;
                logcount = 1;
            } else{
                logcount++;
            }

            posCounter += 7;
        }//end if save D_ij

        if ((x - x2) * (x2 - x1) >= 0.0) //check if end is reached
        {   
            //std::cout << "I'm in this if thingy" << endl;
            for (Int i = 0; i < nvar; i++)
                ystart[i] = y[i];
            if (out.kmax > 0 && abs(out.xsave[out.count - 1] - x2) > 100.0 * abs(x2) * EPS)
                out.save(x, y);
            
            // set particle end position and velocity before exiting.
            //particle.pos[0] = y[0]*GCT::mtopc; particle.pos[1] = y[1]*GCT::mtopc; particle.pos[2] = y[2]*GCT::mtopc;
            //particle.v[0] = y[3]; particle.v[1] = y[4]; particle.v[2] = y[5];

            return;
        }


        if (abs(s.hnext) <= hmin){
            std::cout << "Step size too small in Odeint! \n";
            return;
        }
        h = s.hnext;
    }
    //file.close();
    /*
    std::cout << "Too many steps! \n";
    throw("Too many steps in routine Odeint");
    */
}

    // For generate_samples.cpp
template <class Stepper>
void Odeint<Stepper>::integrate(    std::vector<double> &D_ij, std::vector<double> &D_ij_time, std::vector<double> &r_vect,
                                    Particle &particle, Bfield &bfield, const double x_max)                     //x = t for this program
{
    int count = 0, logcount = 1;                            // Counters to track when to save D_ij
    int posCounter = 0;

    std::array<double, 3> B, r;

    double t;                                               // Used to test if step should be recorded

        // For testing E* and B*
   // Guiding_Center GC;

    //std::ofstream file;
    //file.open("Data/Test_E_B_eff.dat");
    r = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
    bfield.generate_bfield_at_point(x, B, r);

    y[6] = B[0]; y[7] = B[1]; y[8] = B[2];

    derivs(x, y, dydx);
    

    if (dense)
        out.out(-1, x, y, s, h);
    else
        out.save(x, y);

    
    for (nstp = 0; x <= x_max; nstp++)
    {

        if ((x + h * 1.0001 - x2) * (x2 - x1) > 0.0){
            h = x2 - x;
        }

      // Take a step
        s.step(h, derivs);


        if (s.hdid == h){
            ++nok;
        } else{
            ++nbad;
        }
        if (dense){
            out.out(nstp, x, y, s, s.hdid);
        } else{
            out.save(x, y);
        }

        r = { y[0]*GCT::mtopc, y[1]*GCT::mtopc, y[2]*GCT::mtopc };
        bfield.generate_bfield_at_point(x, B, r);

        y[6] = B[0]; y[7] = B[1]; y[8] = B[2];

     //std::cout << "Timestep h: " << h << std::endl;
        double t = x+h;
        if( t/31557600 - (logcount * pow(10, count)) > __DBL_EPSILON__){
            using GCT::mtopc;
            // Add each particles coordinate to the sum in D_ij
            // D_ij holds the upper half, including the diagonal, of a symmetric matrix
            // stored as a 1D-array to comply with MPIs send and recieve functions.
            D_ij[posCounter + 0] += (y[0]*mtopc) * (y[0]*mtopc);                             //D_11
            D_ij[posCounter + 1] += (y[0]*mtopc) * (y[1]*mtopc);                             //D_12
            D_ij[posCounter + 2] += (y[0]*mtopc) * (y[2]*mtopc);                             //D_13
            D_ij[posCounter + 3] += (y[1]*mtopc) * (y[1]*mtopc);                             //D_22
            D_ij[posCounter + 4] += (y[1]*mtopc) * (y[2]*mtopc);                             //D_23
            D_ij[posCounter + 5] += (y[2]*mtopc) * (y[2]*mtopc);                             //D_33
            D_ij[posCounter + 6] += (y[0]*mtopc) * (y[0]*mtopc) + (y[1]*mtopc) * (y[1]*mtopc) + 
                                    (y[2]*mtopc) * (y[2]*mtopc);                             //D_average

            D_ij_time[posCounter / 7] = x/31557600;

            r_vect[3 * (posCounter / 7) + 0] = y[0]*mtopc;
            r_vect[3 * (posCounter / 7) + 1] = y[1]*mtopc;
            r_vect[3 * (posCounter / 7) + 2] = y[2]*mtopc;

            if(logcount == 9){
                count++;
                logcount = 1;
            } else{
                logcount++;
            }

            posCounter += 7;
        }//end if save D_ij

        if ((x - x2) * (x2 - x1) >= 0.0) //check if end is reached
        {   
            //std::cout << "I'm in this if thingy" << endl;
            for (Int i = 0; i < nvar; i++)
                ystart[i] = y[i];
            if (out.kmax > 0 && abs(out.xsave[out.count - 1] - x2) > 100.0 * abs(x2) * EPS)
                out.save(x, y);
            
            // set particle end position and velocity before exiting.
            particle.pos[0] = y[0]; particle.pos[1] = y[1]; particle.pos[2] = y[2];
            particle.v[0] = y[3]; particle.v[1] = y[4]; particle.v[2] = y[5];

            return;
        }


        if (abs(s.hnext) <= hmin){
            std::cout << "Step size too small in Odeint! \n";
            return;
        }
        h = s.hnext;
    }
    //file.close();
    /*
    std::cout << "Too many steps! \n";
    throw("Too many steps in routine Odeint");
    */
}



    // Original integrate function
template <class Stepper>
void Odeint<Stepper>::integrate()           //x = t for this program
{
    int count = 0, logcount = 1;                            //Counters to track when to save D_ij
    int posCounter = 0;
    const double mtopc = 3.24078e-17;
    
    derivs(x, y, dydx);
    

    if (dense)
        out.out(-1, x, y, s, h);
    else
        out.save(x, y);

    
    for (nstp = 0; nstp <= MAXSTP; nstp++)
    {
        //std::cout << "Am I stuck here? \n";

        if ((x + h * 1.0001 - x2) * (x2 - x1) > 0.0){
            h = x2 - x;
        }

        s.step(h, derivs);

        //std::cout << "V_total: " << sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) << endl;

        if (s.hdid == h)
            ++nok;
        else{
            //std::cout << "I did a dumb \n";
            ++nbad;
        }
        if (dense)
            out.out(nstp, x, y, s, s.hdid);
        else
            out.save(x, y);


        if ((x - x2) * (x2 - x1) >= 0.0) //check if end is reached
        {   
            //std::cout << "I'm in this if thingy" << endl;
            for (Int i = 0; i < nvar; i++)
                ystart[i] = y[i];
            if (out.kmax > 0 && abs(out.xsave[out.count - 1] - x2) > 100.0 * abs(x2) * EPS)
                out.save(x, y);
            return;
        }


        if (abs(s.hnext) <= hmin){
            std::cout << "Step size too small in Odeint! \n";
            return;
        }
        h = s.hnext;
    }
   
    /*
    std::cout << "Too many steps! \n";
    throw("Too many steps in routine Odeint");
    */
}

#endif