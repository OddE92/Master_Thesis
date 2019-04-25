#ifndef ODEINT
#define ODEINT

#include "NR3/nr3.h"
#include "Bfield/class_bfield.h"
#include "Trajectory/class_trajectory.h"
#include "Initializer/initializer.h"

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
            D_ij[posCounter + 6] += (y[0]*mtopc) * (y[0]*mtopc) + (y[1]*mtopc) * (y[1]*mtopc) + \
                                    (y[2]*mtopc) * (y[2]*mtopc); //D_average

            D_ij_time[posCounter / 7] = x/31557600;

            r_vect[3 * (posCounter / 7) + 0] = y[0]*mtopc;
            r_vect[3 * (posCounter / 7) + 1] = y[1]*mtopc;
            r_vect[3 * (posCounter / 7) + 2] = y[2]*mtopc;

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
            throw("Step size too small in Odeint");
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
            throw("Step size too small in Odeint");
        }
        h = s.hnext;
    }
   
    /*
    std::cout << "Too many steps! \n";
    throw("Too many steps in routine Odeint");
    */
}

#endif