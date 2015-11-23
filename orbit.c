#include "orbit.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//integration parameters
double const tstart = 0.0;          //time at input of cluster coordinates [Gyr], usually today, i.e. 0.0

double const mdiff = 1.E-7;         //precission
//double const dt0 = 1.E-5;         //initial time-step [Gyr]
double const RMIN = 1.0E-1;         //Smallest allowed separation between galaxies (effectively a softening)
double const dtmax = 0.025;          //maximum time-step [Gyr]
//double const Rgalmin = 10.0;       //minimum galactocentric radius [pc]
//double const Rgalmax = 1.0e10;    //maximum galactocentric radius [pc]
int const VARIABLE_TIMESTEPS = 0;
int const RK4 = 1; // Use a Runge-Kutta? Alt. is leapfrog.

//currently does nothing
void custom_gsl_error_handler(const char * reason,
                              const char * file,
                              int line,
                              int gsl_errno){
    (void)reason;
    (void)file;
    (void)line;
}


int orbit(int ngals,
          struct Params parameters,
          struct Gal *gal,
          struct Snapshot **output_snapshots){


    int ratio;
    double tpast = parameters.tpast;  // add to orbit.h
    double tfuture = parameters.tfuture; // add to orbit.h
    double dt0 = parameters.dt0; // add to orbit.h
    double dtout = parameters.dtout;
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;

    //get position of cluster at t = -tpast
    double sign, tmax, dtoutt, t;
    int err = 0;
    int RECORD_SNAP, WRITE_SNAP;
    if (tpast < 0.0) {
        sign = -1.0;
        tmax = tpast;
        dtoutt = -1.0*dtout;

        // If we don't want to integrate forward, then save snapshots as we go back
        if (tfuture <= 0.0){
            RECORD_SNAP = 1;
            WRITE_SNAP = parameters.snapshot; // Use users choice
        } else {  // otherwise don't save snapshots during backward integration
            RECORD_SNAP = 0;
            WRITE_SNAP = 0;
        }

        t = tstart;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal,
                      parameters, sign, output_snapshots,
                      RECORD_SNAP, WRITE_SNAP);
    }
    if (tfuture > 0.0) {
    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
        sign = 1.0;
        dtoutt = dtout;
        RECORD_SNAP = 1;  // Always save going forward
        WRITE_SNAP = parameters.snapshot;  // Use users choice
        tmax = tfuture;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal,
                      parameters, sign, output_snapshots,
                      RECORD_SNAP, WRITE_SNAP);
    }

    free(gal);
    return err;
}


/* --------------- extrapolation method --------------- */
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
            int WRITE_SNAP){
    int snapnum = 0;
	double tout, diff, dt = 0.0;
	double xe1[3], ve1[3], difftemp;
        //double rt = 1e+5;
    double rt_temp, r, E;
    int k, n, m;
    int err = 0;
	int ngals = parameters.ngals;
	//initialize timesteps
	tout = *t;		    /* time of next output/insertion */
	dt = sign*dt0;                /* initial time step */
	//integrate galaxies
	do {
		/***********
		 * GALAXY *
		 ***********/
        //advance each particle
	    int count = 0;
        int laststep = 0;
	    do {
		    difftemp = 0.0;
	        diff = 0.0;
            // loop over each galaxy
            for (n=0; n<ngals; n++){
                // If galaxy is fixed in place, or stripped then do not advance it
                if ((gal[n].inplace == 1) ||
                    (gal[n].stripped == 1)){
                    continue;
                }
                    // Advance other particles using a fixed or variable time step
		        for (k=0;k<3;k++) {
                    gal[n].post[k] = gal[n].pos[k];
                    gal[n].velt[k] = gal[n].vel[k];
			        xe1[k] = gal[n].pos[k];
			        ve1[k] = gal[n].vel[k];
		        }
                if (VARIABLE_TIMESTEPS) {
		            err = do_step(dt, xe1, ve1, n, gal, parameters);      /* One full step */
		            err = do_step(0.5*dt, gal[n].post, gal[n].velt, n, gal, parameters);  /* Two half steps */
		            err = do_step(0.5*dt, gal[n].post, gal[n].velt, n, gal, parameters);
		            difftemp = sqrt(pow(xe1[0] - gal[n].post[0],2) +
		                            pow(xe1[1] - gal[n].post[1],2) +
		                            pow(xe1[2] - gal[n].post[2],2));
                    if (difftemp > diff) {
                        diff = difftemp;  // hold highest value to compare with mdiff below
                    }
                } else {
		            err = do_step(dt, gal[n].post, gal[n].velt, n, gal, parameters);      /* One full step */
                }
                // error somewhere in integration
                if (err) {
                    return 1;
                }
            }  // end loop over each galaxy

		    if (!VARIABLE_TIMESTEPS || (diff<=mdiff)) {         /* Is difference below accuracy threshold? */
		        *t+=dt;/* If yes -> continue and double step size */
		        //TODO: update the test particles here
                for (n=0; n<ngals; n++){
		            for (k=0;k<3;k++) {
			            gal[n].pos[k]=gal[n].post[k];
			            gal[n].vel[k]=gal[n].velt[k];
		            }
		        }
                // Tidal stripping
                for  (n=0; n<ngals; n++){
                    if ((gal[n].tidal_trunc == 1) &&
                        (gal[n].stripped == 0)) {// tidal truncation turned on - galaxy intact
                        for (m=0; m<ngals; m++){ // look for all galaxies with dynamic friction turned on
                            if ((m != n) && // not self
                            (gal[m].dyn_fric == 1) &&  // dynamical friction on
                            (gal[m].stripped == 0)){  //not stripped
                                r = sqrt(pow(gal[n].pos[0]-gal[m].pos[0], 2) +
                                         pow(gal[n].pos[1]-gal[m].pos[1], 2) +
                                         pow(gal[n].pos[2]-gal[m].pos[2], 2));
                                rt_temp = calc_rt(r, fmin(gal[n].rt, gal[n].r_halo), gal[m], gal[n]);
                                if (rt_temp < 0.0){ // error has occrued
                                    E = binding_energy(gal[n]);
                                    if (E > 0){  // error due to galaxy being stripped
                                        gal[n].stripped = 1;
                                        break;
                                    } else {  // some other error, exit
                                        return 2;
                                    }
                                }
                                if (sign > 0){  // integrating forward
                                    gal[n].rt = fmin(gal[n].rt, rt_temp);
                                } else if (sign < 0){  // integrating backward
                                    gal[n].rt = fmax(gal[n].rt, rt_temp);
                                }

                                if (gal[n].halo_type == 1){ // Dehnen
                                    gal[n].mhalo = gal[n].minit*pow(gal[n].rt/(gal[n].rt+gal[n].r_halo), 3-gal[n].gamma);
                                } else if (gal[n].halo_type == 2){ // NFW
                                    gal[n].mhalo = gal[n].minit*(log(1+(gal[n].rt)/gal[n].r_halo)-gal[n].rt/
                                                                 (gal[n].r_halo+gal[n].rt))/
                                                                (log(1+gal[n].c_halo) -gal[n].c_halo/
                                                                                        (1+gal[n].c_halo));
                                } else { // Plummer
                                    gal[n].mhalo = gal[n].minit*pow(gal[n].rt, 3)/pow(pow(gal[n].r_halo, 2) +
                                                                                      pow(gal[n].rt, 2), 1.5);
                                }
                            }
                        }
                    }
                }
    		    if (VARIABLE_TIMESTEPS) { //If we are with the threshold double the timestep and continue
                    dt = dt*2.0;
                }
		    } else {
		        dt = dt/2.0;  // if we are outside the threshold halve the timestep and try again
		    }
            // Abort if the timestep ever gets too low
		    if (sign*dt < 0.01*dt0 && !laststep) {
		        printf("Aborted... dt = %lf (>%lf, %lf)\n", dt, dt0, sign);
		        return 3;
		    }
		    count++;
                // round to the end of simulation time
    		if (sign*dt > dtmax) {
    		    dt = sign*dtmax;
    		}

		} while (diff>mdiff);       /* Go through loop once and only repeat if difference is too large */
	    if (sign**t>=sign*(tout)) {
            E = binding_energy(gal[n]);
            if (E > 0.0){
                gal[n].stripped = 1;
            }
            //for (n=0; n<ngals; n++) {
            //    if (gal[n].tidal_trunc == 1) {
            //        printf("rt: %10.5f, m: %10.5f, t: %10.5f, E: %10.5f\n", gal[n].rt, gal[n].mhalo, *t, E);
            //    }
            //}

            // First record the snapshot variable
            if (RECORD_SNAP){
                record_snapshot(parameters, gal, *t, snapnum, output_snapshots);
            }
            // Then write out to disk
            if (WRITE_SNAP){
                write_snapshot(parameters, gal, *t, snapnum);
            }
            snapnum += 1;
            tout+=dtout;              /* increase time of output/next insertion */
        }


    } while (sign**t<sign*(tmax));
    // write final snapshot
    if (parameters.snapshot) write_snapshot(parameters, gal, *t, snapnum);
	return 0;

}


/* ---------- advancement ---------- */
int do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, struct Params parameters) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;

    hh = dt*0.5;
    int err = 0;
    if (RK4) {
    	err = getforce_gals(x, v, acc0, gal_num, gal, parameters);
    	for (k=0;k<3;k++) {                /* first half-step */
    	    xt1[k] = *(x+k)+hh**(v+k);
    	    vt1[k] = *(v+k)+hh**(acc0+k);
    	}

    	err = getforce_gals(&xt1[0], &vt1[0], acc1, gal_num, gal, parameters);
    	for (k=0;k<3;k++) {                /* second half-step */
    	    xt2[k] = *(x+k)+hh*vt1[k];
    	    vt2[k] = *(v+k)+hh**(acc1+k);
    	}

    	err = getforce_gals(&xt2[0], &vt2[0], acc2, gal_num, gal, parameters);
    	for (k=0;k<3;k++) {                /* third half-step with results of second half-step */
    	    xt3[k] = *(x+k)+dt*vt2[k];
    	    vt3[k] = *(v+k)+dt**(acc2+k);
    	}

	err = getforce_gals(&xt3[0], &vt3[0], acc3, gal_num, gal, parameters);
	for (k=0;k<3;k++) {                /* Runge-Kutta formula */
	    *(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
	    *(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
	}
    }
    else {  // modified leapfrog
        // ai
        err = getforce_gals(x, v, acc0, gal_num, gal, parameters);
        // vi+1/2 and xi+1/2
        for (k=0;k<3;k++) {
            vt1[k] = *(v+k)+ *(acc0+k)*hh;
            xt1[k] = *(x+k)+ *(v+k)*hh;
        }
        // ai+1/2
        err = getforce_gals(xt1, vt1, acc1, gal_num, gal, parameters);
        //vi+1 and xi+1
        for (k=0;k<3;k++) {
            vt1[k] = *(v+k)+ *(acc1+k)*dt;
            *(x+k) += 0.5*(*(v+k)+vt1[k])*dt;
            *(v+k) = vt1[k];
        }
    }
    return err;
}


int getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, struct Params parameters){

    int i;
    int err = 0;
    double r;
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double vr = 0.0;

    int ngals = parameters.ngals;
    // get the force from all other galaxies
    for (i=0; i<ngals; i++){
        if (i != gal_num){  // skip itself
            //Hernquist bulge
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
                 pow(*(x+1) - gal[i].pos[1], 2) +
                 pow(*(x+2) - gal[i].pos[2], 2));
            if (r > RMIN){
        	    ax += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
        	        (r+gal[i].b1_LMJ))*(*x - gal[i].pos[0])/r;
        	    ay += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
        	        (r+gal[i].b1_LMJ))*(*(x+1) - gal[i].pos[1])/r;
        	    az += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
        	        (r+gal[i].b1_LMJ))*(*(x+2) - gal[i].pos[2])/r;
            }
            //Miyamato disk
    	    r = sqrt(pow(*x - gal[i].pos[0], 2) +
    	             pow(*(x+1) - gal[i].pos[1], 2) +
    	             pow(gal[i].a2_LMJ + sqrt(pow(*(x+2) - gal[i].pos[2], 2)
    	                                      + pow(gal[i].b2_LMJ, 2)), 2)
    	             );
            if (r > RMIN){
        	    ax += -G*gal[i].M2_LMJ/(r*r*r) * (*x - gal[i].pos[0]);
        	    ay += -G*gal[i].M2_LMJ/(r*r*r) * (*(x+1) - gal[i].pos[1]);
        	    az += -G*gal[i].M2_LMJ/(r*r*r) * (gal[i].a2_LMJ + sqrt(pow(*(x+2) - gal[i].pos[2], 2)
        	                                                    + pow(gal[i].b2_LMJ, 2)))
        	                                             / sqrt(pow(*(x+2) - gal[i].pos[2], 2)
        	                                                    + pow(gal[i].b2_LMJ, 2))
        	                                             * (*(x+2) - gal[i].pos[2]);
            }
            // Dark Matter Halo
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
               pow(*(x+1) - gal[i].pos[1], 2) +
               pow(*(x+2) - gal[i].pos[2], 2));
            if (r > RMIN){
                halo_acc(r, gal[i], x, &ax, &ay, &az);
                // dynamical friction
                if (gal[i].dyn_fric == 1) {// is dynamical friction turned on for this galaxy?
                    //relative velocity
                    vx = (*v - gal[i].vel[0]);
                    vy = (*(v+1) - gal[i].vel[1]);
                    vz = (*(v+2) - gal[i].vel[2]);
                    vr = sqrt(vx*vx + vy*vy + vz*vz);

                    err = dynamical_friction(r, vx, vy, vz, vr,
                                             &ax, &ay, &az,
                                             gal[i], gal[gal_num].mhalo,
                                             gal[gal_num].r_halo);

                }
            }
        }
    }
    // update acceleration
    *(a+0) = ax;
    *(a+1) = ay;
    *(a+2) = az;
    return err;
}

// Calculate Dynamical Friction acceleration
int dynamical_friction(double r, double vx, double vy, double vz, double vr,  // orbit velocity and radius
                       double *ax, double *ay, double *az,  // accelerations update in function
                       struct Gal gal,
                       double m_gal, double r_gal){  // companion mass and scale length
    double sigma = 0.0;
    double density = 0.0;
    double dyn_L, dyn_C, dyn_alpha;
    int halo_type = gal.halo_type;
    double mhalo = gal.mhalo;
    double r_halo = gal.r_halo;
    double gamma = gal.gamma;
    double c_halo = gal.c_halo;
    /*
    * COULOMB LOGARITHM
    */
    double coulomb;
    // alternative methods of coulomb logarithm
    if (mhalo/m_gal < 0.2){
        coulomb = 0.0;  // don't have a small galaxy act on a big one
    } else {
        // van den Marel et al. 2012 eq. A1 and discussion in Appendix A
        if (abs(mhalo/m_gal - 1.0) < 0.3){ // within 30% of each other
            dyn_L = dyn_L_eq;
            dyn_C = dyn_C_eq;
            dyn_alpha = dyn_alpha_eq;
        } else {
            dyn_L = dyn_L_uneq;
            dyn_C = dyn_C_uneq;
            dyn_alpha = dyn_alpha_uneq;
        }
        coulomb = fmax(dyn_L, pow(log(r/(dyn_C*r_gal)),
                                    dyn_alpha));
    }
    /*
    * XXXXXXXXXXXXXXXXXXXXXXXXX
    * X and Velocity dispersion
    * XXXXXXXXXXXXXXXXXXXXXXXXX
    */
    halo_sigma(r, gal, &sigma);
    //calculate velocity dispersion

    double X = vr/(sqrt(2.0)*sigma);

    halo_density(r, gal, &density);

    *ax += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vx/pow(vr, 3);

    *ay += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vy/pow(vr, 3);

    *az += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vz/pow(vr, 3);

    return 0;
}

double tidal_condition(double x, void * params)
{
    double *p = (double *)params;
    double MD, MG;
  // p[0] = MD, p[1] = MG, p[2] = aD, p[3] = aG, p[4] = r, p[5] = gammaD, p[6] = gammaG, p[7] = dwarf_type, p[8] = gal_type, p[9] = dwarf_c, p[10] = gal_c
    if (p[7] == 1) { // Dehnen
        MD = p[0]*pow(x/(x+p[2]), 3-p[5]);
    } else if (p[7] == 2){ // NFW
        MD = p[0]*(log(1+(x)/p[2])-(x)/(p[2]+x))/(log(1+p[9]) -p[9]/(1+p[9]));
    } else {// Plummer
        MD = p[0]*pow(x, 3)/pow(p[2]*p[2] + x*x, 1.5);
    }
    if (p[8] == 1){ // Dehnen
        MG = p[1]*pow((p[4]-x)/(p[4]-x+p[3]), 3-p[6]);
    } else if (p[8] == 2){ // NFW
        MG = p[1]*(log(1+(p[4]-x)/p[3])-(p[4]-x)/(p[3]+p[4]-x))/(log(1+p[10]) -p[10]/(1+p[10]));
    }
    else {  //Plummer
        MG = p[1]*pow(p[4]-x, 3)/pow(p[3]*p[3] + pow(p[4]-x, 2), 1.5);
    }

    return fabs(MG/pow(p[4]-x, 3)-
                MD/pow(x, 3));
}

double calc_rt(double r, double rt, struct Gal galG, struct Gal galD)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = rt;
  double a = 1e-1, b = r-1e-1;
  gsl_function F;

  double p[11] = {galD.minit, galG.mhalo, galD.r_halo, galG.r_halo, r,
                  galD.gamma, galG.gamma, galD.halo_type, galG.halo_type,
                  galD.c_halo, galG.c_halo};
  F.function = &tidal_condition;
  F.params = (void *)p;
  gsl_set_error_handler(&custom_gsl_error_handler);

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  status = gsl_min_fminimizer_set (s, &F, m, a, b);
  if (status){
      return -1.0; // GSL reports non-zero error codes
  }

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      if (status) {
          break;
      }

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (a, b, 0.001, 0.001);

    }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status){
      return -1.0; // GSL reports non-zero error codes
  }
  gsl_min_fminimizer_free (s);

  return m;
}


double binding_w(double x, void *params){

    double *p = (double *)params;
    // p[0] = M, p[1] = a, p[2] = gamma, p[3] = c, p[4] = type
    double mass, density;

    if (p[4] == 1) { // Dehnen

        mass = p[0]*pow(x/(x+p[1]), 3-p[2]);
        density = (3.0 - p[2])*p[0]/(4.0*Pi)*
              p[1]/(pow(x, p[2])*pow(x + p[1], 4-p[2]));

    } else if (p[4] == 2){ // NFW

        mass = p[0]*(log(1+(x)/p[1])-(x)/(p[1]+x))/(log(1+p[3]) -p[3]/(1+p[3]));
        density = p[0]/(4.0*Pi*pow(p[1], 3))/
                        (log(1+p[3])-p[3]/(1+p[3]))/
                        (x/p[1]*pow(1+x/p[1], 2));

    } else {// Plummer

        mass = p[0]*pow(x, 3)/pow(p[1]*p[1] + x*x, 1.5);
        density = 3*p[0]/(4*Pi*pow(p[1], 3)*
                pow(1+pow(x/p[1], 2), 2.5));

    }

    return x*density*mass;
}


double binding_t(double x, void *params){

    double *p = (double *)params;
    // p[0] = M, p[1] = a, p[2] = gamma, p[3] = c, p[4] = type
    double density, sigma;

    if (p[4] == 1) { // Dehnen
        sigma = G*p[0]*x*pow(p[1] + x, 3)*
            (-(25.0*pow(p[1], 3) + 52.0*pow(p[1], 2)*x +
                42.0*p[1]*pow(x, 2) + 12.0*pow(x, 3))/
                (12.0*pow(p[1], 4)*pow(p[1] + x, 4)) +
                log((p[1] + x)/x)/pow(p[1], 5));

        density = (3.0 - p[2])*p[0]/(4.0*Pi)*
              p[1]/(pow(x, p[2])*pow(x + p[1], 4-p[2]));

    } else if (p[4] == 2){ // NFW
        // Numerical fit where max is at r=2.16258*a
        double rvmax = 2.16258;
        double VMAX = G*p[0]/
                        (log(1+p[3])-p[3]/(1+p[3]))
                        /(rvmax*p[1])*
                        (log(1+rvmax)-rvmax/(1+rvmax));
        // fitting formula from Zentner and Bullock 2003, eq. 6)
        sigma = 3.0* VMAX *1.4393*pow(x, 0.354)/(1+1.1756*pow(x, 0.725));
        density = p[0]/(4.0*Pi*pow(p[1], 3))/
                        (log(1+p[3])-p[3]/(1+p[3]))/
                        (x/p[1]*pow(1+x/p[1], 2));

    } else {// Plummer
        sigma = pow(p[1], 5)*G*p[0]*
                    pow(1+pow(x/p[1], 2), 2.5)/
                    (6.0*pow(pow(p[1], 2)+pow(x, 2), 3));
        density = 3*p[0]/(4*Pi*pow(p[1], 3)*
                pow(1+pow(x/p[1], 2), 2.5));

    }

    return x*x*density*sigma;
}


double binding_energy(struct Gal gal){
    // init
    int WORKSIZE = 100000;
    double W, T;
    double p[5] = {gal.minit, gal.r_halo, gal.gamma, gal.c_halo, gal.halo_type};
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.params = (void *)p;

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



void write_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber){
    int n;
    int ngals = parameters.ngals;
    char *folder = parameters.outputdir;
    FILE *snapfile;
    char snapname[50];
    double acc0[3];
    sprintf(snapname, "%ssnapshot.csv.%03d", folder, snapnumber);
    snapfile = fopen(snapname, "w");
    fprintf(snapfile,"NAME,X,Y,Z,VX,VY,VZ,AX,AY,AZ,T\n");
    for (n=0; n<ngals; n++){
        getforce_gals(gal[n].pos, gal[n].vel, acc0, n, gal, parameters);
        fprintf(snapfile,"%s,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f\n",
                gal[n].name,
                gal[n].pos[0],
                gal[n].pos[1],
                gal[n].pos[2],
                gal[n].vel[0],
                gal[n].vel[1],
                gal[n].vel[2],
                acc0[0],
                acc0[1],
                acc0[2],
                t);
    }
    fclose(snapfile);
}


void record_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber, struct Snapshot **output_snapshot){
    int n, i;
    int ngals = parameters.ngals;
    for (n=0; n<ngals; n++){
        output_snapshot[snapnumber][n].name = gal[n].name;
        output_snapshot[snapnumber][n].stripped = gal[n].stripped;
        for (i=0; i<3; i++){
            output_snapshot[snapnumber][n].pos[i] = gal[n].pos[i];
            output_snapshot[snapnumber][n].vel[i] = gal[n].vel[i];
        }
        output_snapshot[snapnumber][n].t = t;
    }
}


