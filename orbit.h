#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//constants
#define Pi 3.14159265
#define PI 3.14159265
//#define G  0.0043009211           //gravitational constant in [km^2/s^2/Msun*pc]
//#define G 0.0043021135   // From astropy.constants and units 
#define G 43007.1 // GADGET UNITS!
//#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
//#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define SMALL 1.0E-3
#define SUPERSMALL -1.E50

struct Gal // companion galaxies
{
    double pos[3];
    double vel[3];
    double post[3]; // temp positions during timestep check
    double velt[3]; // temp velocities during timestep check
    int ID;
    double mhalo;
    double minit;
    double r_halo;
    double gamma;
    double a2_LMJ;
    double b2_LMJ;
    double M2_LMJ;
    double M1_LMJ;
    double b1_LMJ;
    double c_halo;
    int dyn_fric;  // is dynamical friction turned on for this galaxy?
    double dyn_C_eq;  // Equal mass terms
    double dyn_L_eq;
    double dyn_alpha_eq;
    double dyn_C_uneq;  // Unequal mass terms
    double dyn_L_uneq;
    double dyn_alpha_uneq;
    int tidal_trunc; // does the galaxy become tidally truncated?
    double rt;
    int halo_type;
    int inplace;
    char *name;
};

struct OrbitStats // orbital statistcs
{
    int apocenters;
    int pericenters;
    int dir;
};

struct Params // Galactic and orbital parameters
{
    double b1_LMJ;        //[pc]
    double M1_LMJ;       //[solar masses]
    double a2_LMJ;       //[pc]
    double b2_LMJ;        //[pc]
    double M2_LMJ;       //[solar masses]
    double Mhalo; //M200 of MW
    double q_halo;  // flattening of halo
    double r_halo; // scale radius of halo
    double gamma; //inner slope of halo
    double c_halo; //NFW concentration
    double halo_type; //NFW = 1, dehnen = 0
    double dyn_C_eq;
    double dyn_L_eq;
    double dyn_alpha_eq;
    double dyn_C_uneq;
    double dyn_L_uneq;
    double dyn_alpha_uneq;
    double tpast;
    double tfuture;
    double dt0;  //snapshot freq
    double dtout;
    int ngals; //number of dwarfs
    char *outputdir; //outputfolder
};

//functions

int orbit(int int_mode,
          int ngals,
          struct Params parameters,
          struct Gal *gal,
          double* output_pos,
          double* output_vel);

int rk4_drv(double *t,
            double tmax,
            double dtout,
            double dt0,
            double mdiff,
            struct Gal *gal,
            struct Params parameters,
            double vorz);

void getforce(double *x, double *v, double *a, struct Params parameters, struct Gal gal);
void getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, struct Params parameters);
void do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, struct Params parameters);

void dynamical_friction(double r, double vx, double vy, double vz, double vr,  // orbit velocity and radius
                        double *ax, double *ay, double *az,  // accelerations update in function
                        int halo_type, double mhalo, double r_halo, double gamma, double c_halo, // Halo properties
                        double dyn_L_eq, double dyn_C_eq, double dyn_alpha_eq, // Roughly equal mass dynamical friction
                        double dyn_L_uneq, double dyn_C_uneq, double dyn_alpha_uneq,  // Unequal mass dynamical friction
                        double m_gal, double r_gal);  // companion mass and friction 

void write_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber);

double tidal_condition (double x, void * params);
double calc_rt(double r, // distance
               double rt, // tidal radius
               struct Gal galG, // galaxy doing the tidal stripping
               struct Gal galD); // galaxy being stripped

//integration parameters
//double const dtout = 5.0;          //time step for output [Myr]
double const tstart = 0.0;          //time at input of cluster coordinates [Myr], usually today, i.e. 0.0
//double const tfuture = 0.0;         //time at end of integration [Myr]
//double const tpast = -6000.0;      //time at beginning of integration [Myr]
double const mdiff = 1.E-7;         //precission
//double const dt0 = 1.E-5;			//initial time-step [Myr]
double const dtmax = 0.025;          //maximum time-step [Myr]
//double const Rgalmin = 10.0;       //minimum galactocentric radius [pc]
//double const Rgalmax = 1.0e10;    //maximum galactocentric radius [pc]
int const VARIABLE_TIMESTEPS = 0;
int const RK4 = 1; // Use a Runge-Kutta? Alt. is leapfrog.

int const tails = 1;                //integrate tidal tail test particles (0= no, 1= yes);
double const Rstop = 20.0;          //increase redge if test particles gets inside r < Rstop, set 0 for no redge parameterscan, else e.g. 20 pc
int const radio = 0;                //say what you're doing
int const tailtype = 0;             //0 = normal, 1 = VL2 maessig, 2 = VL2 besser
double const rtidemax = 1.e9;      //maximum value for rtide

//Write out snapshot during integration?
int const snapshot = 1;
//Compute orbital statistics during integration?
int const orbit_stats = 0;
// which galaxies to test:
int const ref_gal = 2;  // MW, move this to an input parameter later
int const test_gal = 1;  // should be LMC, same as above, move to input

int const DYNAMICALFRICTION_MAIN = 1;  // Compute dynamical friction from fixed galaxy?

