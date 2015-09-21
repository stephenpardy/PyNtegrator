#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

//constants
#define Pi 3.14159265
#define PI 3.14159265
//#define G  0.0043009211           //gravitational constant in [km^2/s^2/Msun*pc]
#define G 43007.1  // GADGET UNITS!
#define SMALL 1.0E-5
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

struct Snapshot  //snapshots to save
{
    char *name;
    double pos[3];
    double vel[3];
    double t;
};

struct Params // Galactic and orbital parameters
{
    double tpast;
    double tfuture;
    double dt0;  //snapshot freq
    double dtout;
    int ngals; //number of dwarfs
    char *outputdir; //outputfolder
};

//functions

int orbit(int ngals,
          struct Params parameters,
          struct Gal *gal,
          struct Snapshot ** output_snapshots);

int rk4_drv(double *t,
            double tmax,
            double dtout,
            double dt0,
            double mdiff,
            struct Gal *gal,
            struct Params parameters,
            double sign,
            struct Snapshot **output_snapshots);

int getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, struct Params parameters);
int do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, struct Params parameters);

int dynamical_friction(double r, double vx, double vy, double vz, double vr,  // orbit velocity and radius
                        double *ax, double *ay, double *az,  // accelerations update in function
                        int halo_type, double mhalo, double r_halo, double gamma, double c_halo, // Halo properties
                        double dyn_L_eq, double dyn_C_eq, double dyn_alpha_eq, // Roughly equal mass dynamical friction
                        double dyn_L_uneq, double dyn_C_uneq, double dyn_alpha_uneq,  // Unequal mass dynamical friction
                        double m_gal, double r_gal);  // companion mass and friction

void write_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber);
void record_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber, struct Snapshot **output_snapshot);

double tidal_condition (double x, void * params);
double calc_rt(double r, // distance
               double rt, // tidal radius
               struct Gal galG, // galaxy doing the tidal stripping
               struct Gal galD); // galaxy being stripped

//integration parameters
double const tstart = 0.0;          //time at input of cluster coordinates [Gyr], usually today, i.e. 0.0

double const mdiff = 1.E-7;         //precission
//double const dt0 = 1.E-5;			//initial time-step [Gyr]
double const dtmax = 0.025;          //maximum time-step [Gyr]
//double const Rgalmin = 10.0;       //minimum galactocentric radius [pc]
//double const Rgalmax = 1.0e10;    //maximum galactocentric radius [pc]
int const VARIABLE_TIMESTEPS = 0;
int const RK4 = 1; // Use a Runge-Kutta? Alt. is leapfrog.

//Write out snapshot during integration?
int const SNAPSHOT = 1;


