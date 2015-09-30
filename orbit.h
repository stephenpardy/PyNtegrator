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
    int stripped;
    double rt;
    int halo_type;
    int inplace;
    char *name;
};

struct Snapshot  //snapshots to save
{
    char *name;
    int stripped;
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
    int snapshot; // save snapshots to disk
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
            struct Snapshot **output_snapshots,
            int RECORD_SNAP,
            int WRITE_SNAP);

int getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, struct Params parameters);
int do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, struct Params parameters);

int dynamical_friction(double r, double vx, double vy, double vz, double vr,  // orbit velocity and radius
                        double *ax, double *ay, double *az,  // accelerations update in function
                        struct Gal gal,
                        double m_gal, double r_gal);  // companion mass and friction
void halo_acc(double r, struct Gal gal, double *x, double *ax, double *ay, double *az);
void halo_sigma(double r, struct Gal gal, double *sigma);
void halo_density(double r, struct Gal gal, double *density);

void write_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber);
void record_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber, struct Snapshot **output_snapshot);

double tidal_condition (double x, void * params);
double calc_rt(double r, // distance
               double rt, // tidal radius
               struct Gal galG, // galaxy doing the tidal stripping
               struct Gal galD); // galaxy being stripped

double binding_energy(struct Gal gal);
double binding_t(double x, void *param);
double binding_w(double x, void *param);
