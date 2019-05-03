#ifndef ODEINT
#define ODEINT

#include "NR3/nr3.h"
#include "Bfield/class_bfield.h"
#include "Trajectory/class_trajectory.h"
#include "Initializer/initializer.h"
#include "Functions/functions.h"

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
        if (!dense)
            throw("dense output not set in Output!");
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
                    Bfield &bfield, Particle &particle, const double x_max);
    void integrate_GC(    Trajectory &trajectory,
                    Bfield &bfield, Particle &particle, const double x_max);
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

template <class Stepper>
void Odeint<Stepper>::integrate_GC(    Trajectory &trajectory,
                                        Bfield &bfield, Particle &particle, const double x_max)         //x = t for this program
{
    int count = 0, logcount = 1;                            // Counters to track when to save D_ij
    int posCounter = 0;

    double t;                                               // Used to test if step should be recorded
    double theta;

    const double mtopc = 3.24078e-17;


    std::vector<double> r(3,0), E_eff_cross_B_hat(3,0), b_hat_prev(3,0), axis(3,0), v_perp_hat(3,0);

    r = trajectory.GC_position;                             // GC_position is initialized in parsec

    bfield.generate_bfield_at_point(0.0, r);

    bfield.calculate_partial_b_hat(bfield.B, trajectory.GC_velocity, trajectory.GC_position, h, x);

    bfield.calculate_B_effective(trajectory.GC_velocity, trajectory.GC_position, trajectory.u);
    bfield.calculate_E_effective(particle, trajectory.v_perp);

    E_eff_cross_B_hat = GCT::vector_cross_product(bfield.E_effective, bfield.B_hat);

    double B_eff_dot_E_eff = GCT::vector_dot_product(bfield.B_effective, bfield.E_effective);

    double B_eff_parallell = GCT::vector_dot_product(bfield.B_hat, bfield.B_effective);

    double B_eff_amp = GCT::vector_amplitude(bfield.B_effective);


        // GC_position is initialized to [pc], but we need [m] during the calculation
 y[0] = trajectory.GC_position[0]/mtopc; y[1] = trajectory.GC_position[1]/mtopc; y[2] = trajectory.GC_position[2]/mtopc;      // y[0-2]
 y[3] = trajectory.u;                                                                                            // y[3]
 y[4] = trajectory.gyrophase;                                                                                    // y[4]
 y[5] = bfield.B_effective[0]; y[6] = bfield.B_effective[1]; y[7] = bfield.B_effective[2];                       // y[5-7]
 y[8] = E_eff_cross_B_hat[0]; y[9] = E_eff_cross_B_hat[1]; y[10] = E_eff_cross_B_hat[2];                         // y[8-9]
 y[11] = B_eff_dot_E_eff;                                                                                        // y[10]
 y[12] = B_eff_parallell;                                                                                        // y[11]
 y[13] = B_eff_amp;                                                                                              // y[12]

std::cout << "y[0-2] at start: " << y[0] << ' ' << y[1] << ' ' << y[2] << std::endl;
    derivs(x, y, dydx);                                     // First step
    

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


        if (s.hdid == h)
            ++nok;
        else{
            ++nbad;
        }
        if (dense)
            out.out(nstp, x, y, s, s.hdid);
        else{
            out.save(x, y);
        }

        b_hat_prev = bfield.B_hat;

        r = { y[0]*mtopc, y[1]*mtopc, y[2]*mtopc };
        bfield.generate_bfield_at_point(x, r);
        
        trajectory.GC_velocity = { y[3]*bfield.B_hat[0], y[3]*bfield.B_hat[1], y[3]*bfield.B_hat[2] };

        //calculate particle.v
        axis = GCT::calculate_vector_rotation(b_hat_prev, bfield.B_hat, theta);     // Axis now holds the normal vector for the rotation
        GCT::rotate_vector_in_plane(trajectory.hat_1, axis, theta);                 // \hat{1} is now rotated to the new plane
        trajectory.a_hat = trajectory.hat_1;                                        // set phase to zero
    //std::cout << "hat_1: " << trajectory.hat_1;
        GCT::rotate_vector_in_plane(trajectory.a_hat, bfield.B_hat, y[4]);          // rotate by the gyrophase
    //std::cout << "a_hat: " << trajectory.a_hat << std::endl;
        v_perp_hat = GCT::vector_cross_product(trajectory.a_hat, bfield.B_hat);
    //std::cout << "v_perp_hat: " << v_perp_hat << std::endl;
        GCT::normalize_vector(v_perp_hat);
        particle.v = {                                                              // Add velocity in perp and parallell directions
            trajectory.v_perp*v_perp_hat[0] + trajectory.v_parallell*bfield.B_hat[0],
            trajectory.v_perp*v_perp_hat[1] + trajectory.v_parallell*bfield.B_hat[1],
            trajectory.v_perp*v_perp_hat[2] + trajectory.v_parallell*bfield.B_hat[2],
        };
/*
        std::cout << "GC_pos: " << y[0] << ' ' << y[1] << ' ' << y[2] << std::endl;
        std::cout << "Particl.v: " << particle.v << std::endl;
        std::cout << "gyrophase: " << y[4] << std::endl;
        std::cout << "u: " << y[3] << std::endl << std::endl;
*/
        bfield.calculate_partial_b_hat(bfield.B, trajectory.GC_velocity, trajectory.GC_position, h, x);

        bfield.calculate_B_effective(trajectory.GC_velocity, r, y[3]);

        bfield.calculate_E_effective(particle, trajectory.v_perp);                              // Assume v_perp = const

        E_eff_cross_B_hat = GCT::vector_cross_product(bfield.E_effective, bfield.B_hat);

        B_eff_dot_E_eff = GCT::vector_dot_product(bfield.B_effective, bfield.E_effective);

        B_eff_parallell = GCT::vector_dot_product(bfield.B_hat, bfield.B_effective);

        B_eff_amp = GCT::vector_amplitude(bfield.B_effective);

        y[5] = bfield.B_effective[0]; y[6] = bfield.B_effective[1]; y[7] = bfield.B_effective[2];
        y[8] = E_eff_cross_B_hat[0]; y[9] = E_eff_cross_B_hat[1]; y[10] = E_eff_cross_B_hat[2];
        y[11] = B_eff_dot_E_eff;
        y[12] = B_eff_parallell;
        y[13] = B_eff_amp;


        double t = x+h;
        if( t/31557600 - (logcount * pow(10, count)) > __DBL_EPSILON__){

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
            std::cout << "I'm in this if thingy" << std::endl;
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

template <class Stepper>
void Odeint<Stepper>::integrate(    std::vector<double> &D_ij, std::vector<double> &D_ij_time, std::vector<double> &r_vect,
                                    Bfield &bfield, Particle &particle, const double x_max)                     //x = t for this program
{
    int count = 0, logcount = 1;                            // Counters to track when to save D_ij
    int posCounter = 0;

    double t;                                               // Used to test if step should be recorded

    const double mtopc = 3.24078e-17;

    std::vector<double> B(3,0), r(3,0);                          // To temporary hold B and r

    // y[6-8] holds the B-field
    // y[0-2] holds the position
    
    r = particle.pos;
    r = { r[0]*mtopc, r[1]*mtopc, r[2]*mtopc };

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

      // Generate B-field at new point
        r = { y[0]*mtopc, y[1]*mtopc, y[2]*mtopc };

        bfield.generate_bfield_at_point(x, B, r);
        y[6] = B[0]; y[7] = B[1]; y[8] = B[2];



        double t = x+h;
        if( t/31557600 - (logcount * pow(10, count)) > __DBL_EPSILON__){

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