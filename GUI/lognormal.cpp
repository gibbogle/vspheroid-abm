#include "mainwindow.h"
#include "log.h"
#include <math.h>
#include <random>
#include <iostream>

LOG_USE();

std::uniform_real_distribution<double> dist(0.0, 1.0);
std::mt19937 mt(12345);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
double random_uni()
{
//    std::random_device rd;
//    std::mt19937 mt(rd());
//    std::uniform_real_distribution<double> dist(0.0, 1.0);

//    std::cout << dist(mt) << "\n";
//    for (int i=0; i<16; ++i)
//        std::cout << dist(mt) << "\n";
//    }

    return dist(mt);
//    return rand()/(RAND_MAX + 1.0);
}

//-----------------------------------------------------------------------------------------
// from: http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
// Implements the Gauss error function.
//   erf(z) = 2 / sqrt(pi) * integral(exp(-t*t), t = 0..z)
//
// fractional error in math formula less than 1.2 * 10 ^ -7.
// although subject to catastrophic cancellation when z in very close to 0
// from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
//-----------------------------------------------------------------------------------------
double erf(double z)
{
   double t = 1.0 / (1.0 + 0.5 * fabs(z));
   // use Horner's method
   double ans = 1 - t * exp( -z*z -  1.26551223 +
         t * ( 1.00002368 +
         t * ( 0.37409196 +
         t * ( 0.09678418 +
         t * (-0.18628806 +
         t * ( 0.27886807 +
         t * (-1.13520398 +
         t * ( 1.48851587 +
         t * (-0.82215223 +
         t * ( 0.17087277))))))))));
   if (z >= 0.0)
       return ans;
   else
       return -ans;
}

//-----------------------------------------------------------------------------------------
// When X is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double pnorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = (x1-mu)/sig;
    e1 = erf(z1/sqrt(2.0))/2;
    z2 = (x2-mu)/sig;
    e2 = erf(z2/sqrt(2.0))/2;
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// When log(X) is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double plognorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = 0;
    z2 = 0;
    if (x1 == 0)
        e1 = -0.5;
    else {
        z1 = (log(x1)-mu)/sig;
        e1 = erf(z1/sqrt(2.0))/2;
    }
    if (x2 == 0)
        e2 = -0.5;
    else {
        z2 = (log(x2)-mu)/sig;
        e2 = erf(z2/sqrt(2.0))/2;
    }
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// Create the lognormal distribution with median = p1, shape = p2
// at n points stored in x[], probability values stored in prob[].
// Note that x[0] = 0.
// The range of x is currently just less than 4*median.  This should be
// OK for values of shape < 2.
// Convert probability into probability density
//-----------------------------------------------------------------------------------------
void MainWindow::create_lognorm_dist(double p1, double p2, int n, double *x, double *prob)
{
    bool use_exp = false;

    if (!use_exp) {
    double xmax, dx, mu_l, sig_l, x1, x2;

    if (p1 >= 0.5)
        xmax = p1*4;
    else
        xmax = p1*8;

    dx = xmax/n;
    mu_l = log(p1);
    sig_l = log(p2);
    for (int ix=0; ix<n; ix++) {
        x1 = (ix - 0.5)*dx;
        x2 = x1 + dx;
        x1 = max(x1,0.0);
        x[ix] = (x1+x2)/2;
        prob[ix] = plognorm(x1,x2,mu_l,sig_l)/(x2-x1);
    }
    } else {
    int nexp = 3;
    double lambda[3];
    double tbase = 20;
    lambda[0] = 1.1;
    lambda[1] = 1.05;
    lambda[2] = 1.1;
    create_expon_dist(tbase,nexp,lambda, n, x, prob);
    }
}

//------------------------------------------------------------------------------------------------------
// Distribution of sum of nexp exponentially distributed RVs, with rate constants lambda[].
// Note: lambda = 1/mean
//------------------------------------------------------------------------------------------------------
void MainWindow::create_expon_dist(double tbase, int nexp, double lambda[], int n, double *x, double *prob)
{
    int i, j, k;
    int nv = 100000;
    double a1, a2, a3, z, dR, U, R;
    if (nexp == 2) {
        dR = (1/lambda[0] + 1/lambda[1])/40;
    } else if (nexp == 3) {
        dR = (1/lambda[0] + 1/lambda[1] + 1/lambda[2])/60;
    }
    for (j=0; j<n; j++) {
        prob[j] = 0;
        x[j] = tbase + j*dR;
    }
    prob[0] = 0;
    // Analytical
    if (nexp == 2) {
        // two cases
        if (lambda[0] == lambda[1]) {
            a1 = lambda[0];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a1*z*exp(-a1*z);
            }
        } else {
            a1 = lambda[0];
            a2 = lambda[1];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a2*(exp(-a1*z) - exp(-a2*z))/(a2-a1);
            }
        }
    } else if (nexp == 3) {
        // three cases: all different, two the same, all the same
        if (lambda[0] == lambda[1] && lambda[0] == lambda[2]) {
            a1 = lambda[0];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a1*a1*z*z*exp(-a1*z)/2;
            }
        } else if (lambda[0] != lambda[1] && lambda[0] != lambda[2] && lambda[1] != lambda[2]) {
            a1 = lambda[0];
            a2 = lambda[1];
            a3 = lambda[2];
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = a1*a2*a3*(exp(-a1*z)/((a2-a1)*(a3-a1)) +
                                    exp(-a2*z)/((a1-a2)*(a3-a2)) +
                                    exp(-a3*z)/((a1-a3)*(a2-a3)));
            }
        } else {
            if (lambda[0] == lambda[1]) {
                a1 = lambda[0];
                a3 = lambda[2];
            } else if (lambda[0] == lambda[2]) {
                a1 = lambda[0];
                a3 = lambda[1];
            } else {
                a1 = lambda[1];
                a3 = lambda[0];
            }
            for (j=1; j<n; j++) {
                z = x[j] - tbase;
                prob[j] = (a1*a1*a3/((a3-a1)*(a3-a1)))*((z*(a3-a1)-1)*exp(-a1*z) + exp(-a3*z));
            }
        }

    }

/*
    // Numerical
    for (j=1; j<n; j++) {
        prob[j] = 0;
        x[j] = tbase + (j+0.5)*dR;
    }
    x[0] = tbase;
    double psum = 0;
    for (k=0; k<nv; k++) {
        R = 0;
        for (i=0; i<nexp; i++) {
            U = random_uni();
            R += (-log(1-U))*lambda[i];
        }
        j = R/dR;
        j = min(j,n-2);
        prob[j+1] += 1;
        psum += 1;
    }
    for (j=1; j<n; j++) {
        prob[j] /= psum;
    }
    prob[0] = 0;
*/
}
