#include "orbit.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


gsl_error_handler_t * gsl_set_error_handler(gsl_error_handler_t * new_handler);


// The acceleration as a function of R
double plummer_acc(double r, struct Gal gal, double x, double x0){
    if ((gal.rt != gal.rt) || (r < gal.rt)){ // normal halo acceleration
        return -G*gal.mhalo / pow(gal.r_halo*gal.r_halo + r*r, 3./2.) * (x - x0);
    } else { // truncated
        return -G*gal.mhalo*(3.*gal.r_halo*gal.r_halo*gal.rt + 2.*pow(gal.rt, 3)) * (x - x0) /
                (gal.r_halo*gal.r_halo * pow(r, 3) * pow(gal.r_halo*gal.r_halo + gal.rt*gal.rt, 3./2.));
    }
}


// Halo velocity dispersion as a function of R
double plummer_sigma(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    return pow(r_halo, 5)*G*mhalo*pow(1.+r*r/(r_halo*r_halo), 5./2.)/(6.*pow(r_halo*r_halo + r*r, 3));
}


// The Halo mass as a function of R
double plummer_mass(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    return mhalo*r*r*r/pow(r_halo*r_halo + r*r, 3./2.);

}


// The halo density as a function of R
double plummer_density(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    return 3.*mhalo/(4.*Pi*r_halo*r_halo*r_halo)*pow(1.+r*r/(r_halo*r_halo), -5./2.);
}


// The W term of the binding energy equation
double plummer_binding_w(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = plummer_density(x, p[0]);
    double mass = plummer_mass(x, p[0]);

    return x*density*mass;
}


// The T term of the binding energy equation
double plummer_binding_t(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = plummer_density(x, p[0]);
    double sigma = plummer_sigma(x, p[0]);

    return x*x*density*sigma;
}


// Compute the binding energy for a galaxy
double plummer_binding_energy(struct Gal gal){
    // init
    gsl_set_error_handler(&custom_gsl_error_handler);

    int WORKSIZE = 100000;
    double W, T;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Integrate for W
    F.function = &plummer_binding_w;
    double result_w, abserr_w;

    int status = gsl_integration_qag(&F, 0, gal.rt, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_w, &abserr_w);
    if (status){
      return 1.0; // GSL reports non-zero error codes
    }
    // Now integrate for T
    F.function = &plummer_binding_t;
    double result_t, abserr_t;

    status = gsl_integration_qag(&F, 0, gal.rt, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_t, &abserr_t);
    if (status){
      return 1.0; // GSL reports non-zero error codes
    }
    //Free integration space
    gsl_integration_workspace_free(workspace);

    W = -4.0*Pi*G*result_w;
    T = 6.0*Pi*result_t;

    return W+T; //dimensionless binding energy - negative is bound
}
