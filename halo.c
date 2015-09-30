#include "orbit.h"
#include <gsl/gsl_math.h>


void halo_acc(double r, struct Gal gal, double *x, double *ax, double *ay, double *az){
    if (gal.halo_type == 1) { // Dehnen
        *ax += -G*gal.mhalo * pow(r/(gal.r_halo + r), -gal.gamma)/
                            pow(gal.r_halo + r, 3) * (*x - gal.pos[0]);
        *ay += -G*gal.mhalo * pow(r/(gal.r_halo + r), -gal.gamma)/
                            pow(gal.r_halo + r, 3) * (*(x+1) - gal.pos[1]);
        *az += -G*gal.mhalo * pow(r/(gal.r_halo + r), -gal.gamma)/
                            pow(gal.r_halo + r, 3) * (*(x+2) - gal.pos[2]);
    } else if (gal.halo_type == 2) {  // NFW
        double constant = -G*gal.mhalo/
                (log(1+gal.c_halo)-gal.c_halo/(1+gal.c_halo));
        *ax += constant * (log(1.0 + r/gal.r_halo)/r -
                        1.0/(gal.r_halo+r)) * (*x - gal.pos[0])/pow(r, 2);
        *ay += constant * (log(1.0 + r/gal.r_halo)/r -
                        1.0/(gal.r_halo+r)) * (*(x+1) - gal.pos[1])/pow(r, 2);
        *az += constant * (log(1.0 + r/gal.r_halo)/r -
                        1.0/(gal.r_halo+r)) * (*(x+2) - gal.pos[2])/pow(r, 2);

    } else { // Plummer sphere
        *ax += -2.0*G*gal.mhalo* (*x - gal.pos[0])/
            pow(pow(gal.r_halo, 2)+pow(r, 2), 1.5);
        *ay += -2.0*G*gal.mhalo* (*(x+1) - gal.pos[1])/
            pow(pow(gal.r_halo, 2)+pow(r, 2), 1.5);
        *az += -2.0*G*gal.mhalo* (*(x+2) - gal.pos[2])/
            pow(pow(gal.r_halo, 2)+pow(r, 2), 1.5);
    }
}


void halo_sigma(double r, struct Gal gal, double *sigma){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double c_halo = gal.c_halo;
    int halo_type = gal.halo_type;
    if (halo_type == 1){
    // Hernquist 1D velocity dispersion - CForm from Mathematica -- only good for gamma == 1
        *sigma = 3.0*sqrt(G*mhalo*r*pow(r_halo + r, 3)*
            (-(25.0*pow(r_halo, 3) + 52.0*pow(r_halo, 2)*r +
                42.0*r_halo*pow(r, 2) + 12.0*pow(r, 3))/
                (12.0*pow(r_halo, 4)*pow(r_halo + r, 4)) +
                log((r_halo + r)/r)/pow(r_halo, 5)));
    } else if (halo_type == 2) {  //NFW
        // Numerical fit where max is at r=2.16258*a
        double rvmax = 2.16258;
        double VMAX = sqrt(G*mhalo/
                        (log(1+c_halo)-c_halo/(1+c_halo))
                        /(rvmax*r_halo)*
                        (log(1+rvmax)-rvmax/(1+rvmax)));
        // fitting formula from Zentner and Bullock 2003, eq. 6)
        *sigma = 3.0* VMAX *1.4393*pow(r, 0.354)/(1+1.1756*pow(r, 0.725));
    } else { //Plummer
        *sigma = 3.0*sqrt(pow(r_halo, 5)*G*mhalo*
                    pow(1+pow(r/r_halo, 2), 2.5)/
                    (6.0*pow(pow(r_halo, 2)+pow(r, 2), 3)));
    }
}


void halo_density(double r, struct Gal gal, double *density){
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double c_halo = gal.c_halo;
    double gamma = gal.gamma;
    int halo_type = gal.halo_type;
    if (halo_type == 1) { //Hernquist
        *density = (3.0 - gamma)*mhalo/(4.0*Pi)*
              r_halo/(pow(r, gamma)*pow(r + r_halo, 4-gamma));
    } else if (halo_type == 2) { //NFW
        // where rho0 = Mvir/(4*Pi*a^3)/(log(1+c)-c/(1+c))
        *density = mhalo/(4.0*Pi*pow(r_halo, 3))/
                        (log(1+c_halo)-c_halo/(1+c_halo))/
                        (r/r_halo*pow(1+r/r_halo, 2));
    } else { //Plummer
        *density = 3*mhalo/(4*Pi*pow(r_halo, 3)*
                pow(1+pow(r/r_halo, 2), 2.5));
    }
}