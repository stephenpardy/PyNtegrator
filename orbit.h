#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//constants
#define Pi 3.14159265
#define PI 3.14159265
#define G  0.0043009211           //gravitational constant in [km^2/s^2/Msun*pc]
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define SMALL 1.0E-3
#define SUPERSMALL -1.E50

struct Gal // companion galaxies
{
    double pos[3];
    double vel[3];
    double post[3]; // temp positions during timestep check
    double velt[3]; // temp velocities durin timestep check
    int ID;
    double mhalo;
    double r_halo;
};

//functions

int orbit(PyDictObject *parameters);

int rk4_drv(double *t,
            double tmax,
            double dtout,
            double mdiff,
            struct Gal *gal,
            double vorz);

void getforce(double *x, double *v, double *a);
void getforce_gals(double *x, double *v, double *a, struct Gal *gal);
void do_step(double dt, double *x, double *v, struct Gal *gal);
double get_gauss(void);
void convert(double *x,
             double *v,
             double *dsun,
             double *vrsun,
             double *vr,
             double *l,
             double *b,
             double *lcosb,
             double *RA,
             double *DEC,
             double *mu_alpha,
             double *mu_alphacosdelta,
             double *mu_delta,
             double *mutemp,
             double *PAtemp,
             int coordtype_coco,
             int vcoordtype,
             int radiococo,
             double vLSRtemp,
             double rgalsun);

void shellsort_reverse_1d(double *array, int N);
void shellsort_1d(double *array, int N);
void shellsort(double **array, int N, int k);
void shellsort_reverse(double **array, int N, int k);
double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);

//integration parameters
double const dtout = 15.0;          //time step for output [Myr]
double const tstart = 0.0;          //time at input of cluster coordinates [Myr], usually today, i.e. 0.0
double const tfuture = 0.0;         //time at end of integration [Myr]
//double const tpast = -6000.0;      //time at beginning of integration [Myr]
double const mdiff = 1.E-4;         //precission
double const dt0 = 1.E-5;			//initial time-step [Myr]
double const dtmax = 50.0;          //maximum time-step [Myr]
double const Rgalmin = 10.0;       //minimum galactocentric radius [pc]
double const Rgalmax = 1.0e10;    //maximum galactocentric radius [pc]

int const tails = 1;                //integrate tidal tail test particles (0= no, 1= yes);
double const Rstop = 20.0;          //increase redge if test particles gets inside r < Rstop, set 0 for no redge parameterscan, else e.g. 20 pc
int const radio = 0;                //say what you're doing
int const tailtype = 0;             //0 = normal, 1 = VL2 maessig, 2 = VL2 besser
double const rtidemax = 1.e9;      //maximum value for rtide

//Write out snapshot after each integration?
int const snapshot = 0;
double const dtsnap = 5.0;           //timestep to save snapshot [Myr]

int const ngals = 2;  // Number of dwarfs

//potential parameters
int const gpot = 3;             //type of Galactic potential (1= Allen & Santillan (1991), 2= log-halo (Koposov et al.), 3= NFW (Irrgang et al.))

//Allen & Santillan potential constants
double const b1 = 230.0;        //I12 //[pc]
double const M1 = 9.4888e09;    //I12 //[solar masses]
double const a2 = 4220.0;       //I12 //[pc]
double const b2 = 292.0;        //I12 //[pc]
double const M2 = 6.62592e10;   //I12 //[solar masses]
double const a3 = 6834.07;      //I12 //[pc]
double const M3 = 2.36176e10;   //I12 //[solar masses]
double const qz = 1.0;          //halo flattening along z axis

//Koposov et al. (2010) potential constants
//nothing to specify; Vc = sqrt(GM3/a3)

//Law, Majewski & Johnston (2009) potential constants
double const b1_LMJ = 700.0;        //[pc]
double const M1_LMJ = 3.4e10;       //[solar masses]
double const a2_LMJ = 6500.0;       //[pc]
double const b2_LMJ = 260.0;        //[pc]
double const M2_LMJ = 1.0e11;       //[solar masses]

//Galactic North Pole parameters
double const alphaGNP = 192.859508; //Galactic north pole in J2000 coordinates
double const deltaGNP = 27.128336;
double const PAGNP = 122.932;       //Position angle with respect to equatorial pole

//Halo parameters
double const Mhalo = 1.5e12; //M200 of MW
double const q_halo = 1.0;  // flattening of halo
double const r_halo = 3e+4; // scale radius of halo... I think

//solar parameters
double const rgalsun = 8330.0;      //solar Galactocentric radius [pc] (standard = 8330.0; Gillessen et al. 2009)
double const vLSR = 239.5;          //rotational velocity of local standard of rest (Reid & Brunthaler 2004 - vysun)
//double vxsun = 0.0;
//double vysun = 0.0;
//double vzsun = 0.0;
//double vxsun = 8.83;              //±0.24 - solar motion with respect to the LSR from Coskunoglu et al. (2011) [km/s]
//double vysun = 14.19;             //±0.34 - 20453 RAVE stars
//double vzsun = 6.57;              //±0.21
double const vxsun = 11.1;          //+0.69/−0.75 - solar motion with respect to the LSR from Schoenrich, Binney & Dehnen (2010) [km/s]
double const vysun = 12.24;         //+0.47−0.47 //24.0;// Bovy et al. (2012)
double const vzsun = 7.25;          //+0.37−0.36
//double vxsun = 10.0;              //solar motion with respect to the LSR from Dehnen & Binney (1998) [km/s]
//double vysun = 5.3;
//double vzsun = 7.2;
//double vxsun = 10.4;              //solar motion with respect to the LSR from Johnson & Soderblom (1987) [km/s]
//double vysun = 14.8;
//double vzsun = 7.3;
//double vxsun = 11.0;              //solar motion with respect to the LSR from Ratnatunga, Bahcall & Casertano (1989) [km/s]
//double vysun = 14.0;
//double vzsun = 7.5;

//fixed cluster parameters
double const vr =  106.7;//110.7;//Harris       //NGC 5466 //heliocentric radial velocity [km/s]
double const redge = 20.0;       //cluster edge radius [pc] CHANGE THIS
double const vrexcess = 0.0;    //mean escape velocity ~ 0.5-2 times velocity dispersion [km/s]
double const rexcess = 0.0;     //displacement from Lagrange point [rtide]
double const R4 = 20.0;         //cluster plummer radius CHANGE THIS
//double const b = 73.59;              //NGC 5466 //galactic lattitude [deg]
//double const l = 42.15;               //NGC 5466 //galactic longitude [deg]  (provide l or l*cos(b), set other to 0.0)
//double const lcosb = 0.0;             //NGC 5466 //galactic longitude times cosine of galactic lattitude [deg]

//double const vr =  -58.7;       //Pal 5 //heliocentric radial velocity [km/s]
//double const redge = 20.0;       //cluster edge radius [pc]
//double const vrexcess = 0.0;    //mean escape velocity ~ 0.5-2 times velocity dispersion [km/s]
//double const rexcess = 0.0;     //displacement from Lagrange point [rtide]
//double const R4 = 20.0;         //cluster plummer radius
//double const b = 45.860;              //Pal5 //galactic lattitude [deg]
//double const l = 0.852;               //Pal5 //galactic longitude [deg]  (provide l or l*cos(b), set other to 0.0)
//double const lcosb = 0.0;             //Pal5 //galactic longitude times cosine of galactic lattitude [deg]
