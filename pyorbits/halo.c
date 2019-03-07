#include "orbit.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


gsl_error_handler_t * gsl_set_error_handler(gsl_error_handler_t * new_handler);


// The acceleration as a function of R
double halo_acc(double r, struct Gal gal, double x, double x0){
    if ((gal.rt != gal.rt) || (r < gal.rt)){ // normal halo acceleration
        return -G*gal.mhalo * pow(r, -gal.gamma)*
                        pow(gal.r_halo + r, -3+gal.gamma) * (x - x0);
    } else { // truncated
        return -G*gal.mhalo * pow(gal.rt/(gal.r_halo + gal.rt), 3 -gal.gamma)*
                        (x - x0)/pow(r, 3);
    }
}


// Sigma term to solve in the halo_sigma function
double sigma_term(double x, void *params){
    struct Gal *p = (struct Gal *)params;
    return pow(x, 1-2*p[0].gamma)/pow(x+p[0].r_halo, 7-2*p[0].gamma);
}


// Halo velocity dispersion as a function of R
double halo_sigma(double r, struct Gal gal){
    // init
    gsl_set_error_handler(&custom_gsl_error_handler);
    int WORKSIZE = 100000;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Now integrate for T
    F.function = &sigma_term;
    double result_s, abserr_s;

    int status = gsl_integration_qag(&F, r, 1e3, 0.0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_s, &abserr_s);
    if (status){
        return -1.0;
    }
    //Free integration space
    gsl_integration_workspace_free(workspace);

    return result_s*G*gal.mhalo*pow(r, gal.gamma)*pow(r + gal.r_halo, 4-gal.gamma);
}


// The sigma term in a turncated halo, used below
double sigma_term_trunc(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = halo_density(x, p[0]);
    double acc = halo_acc(x, p[0], x, 0.0);

    return density*acc;
}


// The halo velocity dispersion as a function of R, in a truncated halo
double halo_sigma_trunc(double r, struct Gal gal){
    // init
    gsl_set_error_handler(&custom_gsl_error_handler);
    int WORKSIZE = 100000;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Now integrate for T
    F.function = &sigma_term_trunc;
    double result_s, abserr_s;

    int status = gsl_integration_qag(&F, r, 1e3, 0.0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_s, &abserr_s);
    if (status){
        return -1.0;
    }
    //Free integration space
    gsl_integration_workspace_free(workspace);

    double density = halo_density(r, gal);

   return result_s/density;
}


// The Halo mass as a function of R
double halo_mass(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double gamma = gal.gamma;
    return mhalo*pow(r/(r+r_halo), 3-gamma);

}


// The halo density as a function of R
double halo_density(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double gamma = gal.gamma;
    return (3.0 - gamma)*mhalo/(4.0*Pi)*
            r_halo/(pow(r, gamma)*pow(r + r_halo, 4-gamma));
}


// The W term of the binding energy equation
double binding_w(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = halo_density(x, p[0]);
    double mass = halo_mass(x, p[0]);

    return x*density*mass;
}


// The T term of the binding energy equation
double binding_t(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = halo_density(x, p[0]);
    double sigma = halo_sigma(x, p[0]);

    return x*x*density*sigma;
}


// Compute the binding energy for a galaxy
double binding_energy(struct Gal gal){
    // init
    gsl_set_error_handler(&custom_gsl_error_handler);

    int WORKSIZE = 100000;
    double W, T;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Integrate for W
    F.function = &binding_w;
    double result_w, abserr_w;

    int status = gsl_integration_qag(&F, 0, gal.rt, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_w, &abserr_w);
    if (status){
      return 1.0; // GSL reports non-zero error codes
    }
    // Now integrate for T
    F.function = &binding_t;
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


// Compute the mass growth of a galaxy
double mass_growth(double t, struct Gal gal){
    //Aquarius mass growth function
    double z = -0.843 * log(1 - (-1*t)/11.32);  // fitting formula
    return gal.minit * pow(1 + z, 2.23) *
                       exp(-4.49*(sqrt(1 + z) - 1.0));  // Aquarius fitting formula for MW
}
