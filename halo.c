#include "orbit.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>

double halo_acc(double r, struct Gal gal, double x, double x0){
    if ((gal.rt != gal.rt) || (r < gal.rt)){ // normal halo acceleration
        return -G*gal.mhalo * pow(r, -gal.gamma)*
                        pow(gal.r_halo + r, -3+gal.gamma) * (x - x0);
    } else { // truncated
        return -G*gal.mhalo * pow(gal.rt/(gal.r_halo + gal.rt), 3 -gal.gamma)*
                        (x - x0)/pow(r, 3);
    }
}


double sigma_term(double x, void *params){
    struct Gal *p = (struct Gal *)params;
    return pow(x, 1-2*p[0].gamma)/pow(x+p[0].r_halo, 7-2*p[0].gamma);
}


double halo_sigma(double r, struct Gal gal){
    // init
    int WORKSIZE = 100000;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Now integrate for T
    F.function = &sigma_term;
    double result_s, abserr_s;

    gsl_integration_qag(&F, r, 1e3, 0.0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_s, &abserr_s);

    //Free integration space
    gsl_integration_workspace_free(workspace);

    return result_s*G*gal.mhalo*pow(r, gal.gamma)*pow(r + gal.r_halo, 4-gal.gamma);
}


double sigma_term_trunc(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = halo_density(x, p[0]);
    double acc = halo_acc(x, p[0], x, 0.0);

    return density*acc;
}


double halo_sigma_trunc(double r, struct Gal gal){
    // init
    int WORKSIZE = 100000;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Now integrate for T
    F.function = &sigma_term_trunc;
    double result_s, abserr_s;

    gsl_integration_qag(&F, r, 1e3, 0.0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_s, &abserr_s);

    //Free integration space
    gsl_integration_workspace_free(workspace);

    double density = halo_density(r, gal);

   return result_s/density;
}


double halo_mass(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double gamma = gal.gamma;
    return mhalo*pow(r/(r+r_halo), 3-gamma);

}


double halo_density(double r, struct Gal gal){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double gamma = gal.gamma;
    return (3.0 - gamma)*mhalo/(4.0*Pi)*
            r_halo/(pow(r, gamma)*pow(r + r_halo, 4-gamma));
}


double binding_w(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = halo_density(x, p[0]);
    double mass = halo_mass(x, p[0]);

    return x*density*mass;
}


double binding_t(double x, void *params){

    struct Gal *p = (struct Gal *)params;
    double density = halo_density(x, p[0]);
    double sigma = halo_sigma(x, p[0]);

    return x*x*density*sigma;
}


double binding_energy(struct Gal gal){
    // init
    int WORKSIZE = 100000;
    double W, T;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = &gal;

    // Integrate for W
    F.function = &binding_w;
    double result_w, abserr_w;

    gsl_integration_qag(&F, 0, gal.rt, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_w, &abserr_w);

    // Now integrate for T
    F.function = &binding_t;
    double result_t, abserr_t;

    gsl_integration_qag(&F, 0, gal.rt, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result_t, &abserr_t);

    //Free integration space
    gsl_integration_workspace_free(workspace);

    W = -4.0*Pi*G*result_w;
    T = 6.0*Pi*result_t;

    return W+T; //dimensionless binding energy - negative is bound
}
