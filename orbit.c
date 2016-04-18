#include "orbit.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>


//integration parameters
double const tstart = 0.0;          //time at input of cluster coordinates [Gyr], usually today, i.e. 0.0

double const mdiff = 1.E-7;         //precission
double const RMIN = 1e-4;         //Smallest allowed separation between galaxies (effectively a softening)
double const dtmax = 0.025;          //maximum time-step [Gyr]

//currently does nothing
void custom_gsl_error_handler(const char * reason,
                              const char * file,
                              int line,
                              int gsl_errno){

    fprintf(stdout, "GSL Error %s Number: %d:\n", reason, gsl_errno);
    fprintf(stdout, "Error occured in file %s at line %d\n", file, line);

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
    int VARIABLE_TIMESTEPS = parameters.variabletimesteps;
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;

    //get position of cluster at t = -tpast
    double sign, tmax, dtoutt, t;
    int err = 0;
    int RECORD_SNAP, WRITE_SNAP, WRITE_TRACERS;
    if (tpast < 0.0) {
        sign = -1.0;
        tmax = tpast;
        dtoutt = -1.0*dtout;

        // If we don't want to integrate forward, then save snapshots as we go back
        if (tfuture <= 0.0){
            RECORD_SNAP = 1;
            WRITE_SNAP = parameters.snapshot; // Use users choice
            WRITE_TRACERS = parameters.write_tracers; // Use users choice
            // Init tracers
            init_tracers(gal, ngals);

        } else {  // otherwise don't save snapshots during backward integration
            RECORD_SNAP = 0;
            WRITE_SNAP = 0;
            WRITE_TRACERS = 0;
        }

        t = tstart;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal,
                      parameters, sign, output_snapshots, VARIABLE_TIMESTEPS,
                      RECORD_SNAP, WRITE_SNAP, WRITE_TRACERS);
        if (err > 0){
            return err;
        }
    }
    if (tfuture > 0.0) {
    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
        sign = 1.0;
        dtoutt = dtout;
        RECORD_SNAP = 1;  // Always save going forward
        WRITE_SNAP = parameters.snapshot;  // Use users choice
        WRITE_TRACERS = parameters.write_tracers;  // Use users choice
        // Init tracers
        init_tracers(gal, ngals);
        tmax = tfuture;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal,
                      parameters, sign, output_snapshots, VARIABLE_TIMESTEPS,
                      RECORD_SNAP, WRITE_SNAP, WRITE_TRACERS);
    }

    //free(gal);
    return err;
}

void init_tracers(struct Gal *gal, int ngals){
srand(100);
for (int n=0; n<ngals; n++){
    if (gal[n].test_particles.nparticles > 0){
        for (int i=0; i<1000; i++){
            for(int j=0; j<3; j++){
                gal[n].test_particles.pos[i*3+j] = gal[n].pos[j]+gal[n].pos[j]*(2.0*(rand() / ((double)RAND_MAX + 1))-1.0)/100.;
                gal[n].test_particles.vel[i*3+j] = gal[n].vel[j];
            }
        }
    }
}
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
            int VARIABLE_TIMESTEPS,
            int RECORD_SNAP,
            int WRITE_SNAP,
            int WRITE_TRACERS){
    int snapnum = 0;
	double tout, diff, dt = 0.0;
	double xe1[3], ve1[3], difftemp;
        //double rt = 1e+5;
    double rt_temp, r, E;
    float z;
    //int k, n, m;
    int err = 0;
	int ngals = parameters.ngals;
	//initialize timesteps
	tout = *t;		    /* time of next output/insertion */
	dt = sign*dt0;                /* initial time step */
	//integrate galaxies
	do {
        //advance each particle
	    int count = 0;
        int laststep = 0;
	    do {
		    difftemp = 0.0;
	        diff = 0.0;
            // loop over each galaxy
            for (int n=0; n<ngals; n++){
                // If galaxy is fixed in place, or stripped then do not advance it
                if ((gal[n].inplace == 1) ||
                    (gal[n].stripped == 1)){
                    continue;
                }
                    // Advance other particles using a fixed or variable time step
		        for (int k=0; k<3; k++) {
                    gal[n].post[k] = gal[n].pos[k];
                    gal[n].velt[k] = gal[n].vel[k];
			        xe1[k] = gal[n].pos[k];
			        ve1[k] = gal[n].vel[k];
		        }
                if (VARIABLE_TIMESTEPS) {
		            err = do_step(dt, xe1, ve1, n, gal, ngals);      /* One full step */
		            err = do_step(0.5*dt, gal[n].post, gal[n].velt, n, gal, ngals);  /* Two half steps */
		            err = do_step(0.5*dt, gal[n].post, gal[n].velt, n, gal, ngals);
		            difftemp = sqrt(pow(xe1[0] - gal[n].post[0],2) +
		                            pow(xe1[1] - gal[n].post[1],2) +
		                            pow(xe1[2] - gal[n].post[2],2));
                    if (difftemp > diff) {
                        diff = difftemp;  // hold highest value to compare with mdiff below
                    }
                } else {
		            err = do_step(dt, gal[n].post, gal[n].velt, n, gal, ngals);      /* One full step */
                }
                // error somewhere in integration
                if (err) {
                    return 1;
                }
            }  // end loop over each galaxy

		    if (!VARIABLE_TIMESTEPS || (diff<=mdiff)) {         /* Is difference below accuracy threshold? */
		        *t+=dt;/* If yes -> continue and double step size */
		        //Update the test particles
                for (int n=0; n<ngals; n++){
                    if (gal[n].test_particles.nparticles > 0){  // if we have test particles
                        err = do_step_tracers(dt, &gal[n].test_particles, gal, ngals);
                        //diagnostic
                        /*
                        for (int m=0; m<gal[n].test_particles.nparticles; m++){
                            printf("ERR: %d (Pos, Vel): ", err);
                            for(int i=0; i<3; i++){
                                printf("(%g, %g) ",
                                        gal[n].test_particles.pos[m*3+i],
                                        gal[n].test_particles.vel[m*3+i]);
                            }
                            printf("\n");
                        }
                        */
                        // end diagnostic
                    }
                }
                // Update galaxy positions
                for (int n=0; n<ngals; n++){
		            for (int k=0;k<3;k++) {
			            gal[n].pos[k]=gal[n].post[k];
			            gal[n].vel[k]=gal[n].velt[k];
		            }
		        }
                //printf("\n");
                // Tidal stripping
                for (int n=0; n<ngals; n++){
                    if ((gal[n].tidal_trunc == 1) &&
                        (gal[n].stripped == 0)) {// tidal truncation turned on - galaxy intact
                        for (int m=0; m<ngals; m++){ // look for all galaxies with dynamic friction turned on
                            if ((m != n) && // not self
                            (gal[m].dyn_fric == 1) &&  // dynamical friction on
                            (gal[m].stripped == 0)){  //not stripped
                                r = sqrt(pow(gal[n].pos[0]-gal[m].pos[0], 2) +
                                         pow(gal[n].pos[1]-gal[m].pos[1], 2) +
                                         pow(gal[n].pos[2]-gal[m].pos[2], 2));
                                rt_temp = calc_rt(r, fmin(gal[n].rt, gal[n].r_halo), gal[m], gal[n]);
                                if (rt_temp < 0.0){ // error has occured
                                    E = binding_energy(gal[n]);
                                    if (E > 0.0){  // error due to galaxy being stripped
                                        gal[n].stripped = 1;
                                        printf("Galaxy %s stripped at time %g \n", gal[n].name, *t);
                                        break;
                                    } else {  // some other error, exit
                                        printf("problem (%g) with: %s (%s), E: %g, r: %g, rt: %g\n", rt_temp,
                                                                                       gal[n].name,
                                                                                       gal[m].name,
                                                                                       E, r, gal[n].rt);
                                        return 2;
                                    }
                                }
                                if (sign > 0){  // integrating forward
                                    gal[n].rt = fmin(gal[n].rt, rt_temp);
                                } else if (sign < 0){  // integrating backward
                                    gal[n].rt = fmax(gal[n].rt, rt_temp);
                                }

                                //gal[n].mtidal = halo_mass(gal[n].rt, gal[n]);
                            }
                        }
                    }
                }
                if (*t < 0){  // only grow when in the past
                    for (int n=0; n<ngals; n++){
                        if (gal[n].mass_growth == 1){
                            //Aquarius mass growth function
                            z = -0.843 * log(1 - (-1**t)/11.32);  // fitting formula
                            gal[n].mhalo = gal[n].minit *
                                            pow(1 + z, 2.23) *
                                            exp(-4.49*(sqrt(1 + z) - 1.0));  // Aquarius fitting formula
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
                printf("Aborted... dt = %g (>%lf, %lf)\n", dt, dt0, sign);
                return 3;
            }
            count++;
            // round to the max allowed timestep
            if (sign*dt > dtmax) {
                dt = sign*dtmax;
            }

            // round to the end of simulation
            if (sign*dt+sign*(*t) > sign*(tmax)) {
                dt = tmax - *t;
            }

		} while (diff>mdiff);       /* Go through loop once and only repeat if difference is too large */
        //for (int n=0; n<ngals; n++) {
        if (sign**t>=sign*(tout)) {
            /*DIAGNOSTIC
            //E = binding_energy(gal[n]);
            //if (E > 0.0){
            //    gal[n].stripped = 1;
            //}
            //for (int n=0; n<ngals; n++) {
            //    if (gal[n].tidal_trunc == 1) {
            //        printf("rt: %10.5f, m: %10.5f, t: %10.5f, E: %10.5f\n", gal[n].rt, gal[n].mtidal, *t, E);
            //    }
            }
            END DIAGNOSTIC*/

            // First record the snapshot variable
            if (RECORD_SNAP){
                record_snapshot(ngals, gal, *t, snapnum, output_snapshots);
            }
            // Then write out to disk
            if (WRITE_SNAP){
                write_snapshot(parameters, gal, *t, snapnum);
            }
            // Then write tracers
            if (WRITE_TRACERS){
                write_tracers(parameters, gal, *t, snapnum);
            }

            snapnum += 1;
            tout+=dtout;              /* increase time of output/next insertion */
        }


    } while (sign**t<sign*(tmax));
    // write final snapshot
    if (WRITE_SNAP) write_snapshot(parameters, gal, *t, snapnum);
	return 0;

}


/* ---------- advancement ---------- */
int do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, int ngals) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;

    hh = dt*0.5;
    int err = 0;
	err = getforce_gals(x, v, acc0, gal_num, gal, ngals);
    if (err) {
        return err;
    }
    for (k=0;k<3;k++) {                /* first half-step */
	    xt1[k] = *(x+k)+hh**(v+k);
	    vt1[k] = *(v+k)+hh**(acc0+k);
	}

	err = getforce_gals(&xt1[0], &vt1[0], acc1, gal_num, gal, ngals);
    if (err) {
        return err;
    }
    for (k=0;k<3;k++) {                /* second half-step */
	    xt2[k] = *(x+k)+hh*vt1[k];
	    vt2[k] = *(v+k)+hh**(acc1+k);
	}

	err = getforce_gals(&xt2[0], &vt2[0], acc2, gal_num, gal, ngals);
    if (err) {
        return err;
    }
	for (k=0;k<3;k++) {                /* third half-step with results of second half-step */
	    xt3[k] = *(x+k)+dt*vt2[k];
	    vt3[k] = *(v+k)+dt**(acc2+k);
	}

	err = getforce_gals(&xt3[0], &vt3[0], acc3, gal_num, gal, ngals);
    if (err) {
        return err;
    }
	for (k=0;k<3;k++) {                /* Runge-Kutta formula */
	    *(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
	    *(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
	}

    return err;
}


/* ---------- advancement of test particles ---------- */
int do_step_tracers(double dt, struct Tracer *test_particles, struct Gal *gal, int ngals) {
    double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];

    hh = dt*0.5;
    int err = 0;
    double x[3];
    double v[3];
    for (int n=0; n<(*test_particles).nparticles; n++){
        for (int k=0; k<3; k++){
            x[k] = (*test_particles).pos[3*n+k];
            v[k] = (*test_particles).vel[3*n+k];
        }

        err = getforce_tracers(x, v, acc0, gal, ngals);
        for (int k=0; k<3; k++) {                /* first half-step */
            xt1[k] = x[k]+hh*v[k];
            vt1[k] = v[k]+hh*acc0[k];
        }

        err = getforce_tracers(&xt1[0], &vt1[0], acc1, gal, ngals);
        for (int k=0; k<3; k++) {                /* second half-step */
            xt2[k] = x[k]+hh*vt1[k];
            vt2[k] = v[k]+hh*acc1[k];
        }

        err = getforce_tracers(&xt2[0], &vt2[0], acc2, gal, ngals);
        for (int k=0; k<3; k++) {                /* third half-step with results of second half-step */
            xt3[k] = *(x+k)+dt*vt2[k];
            vt3[k] = *(v+k)+dt*acc2[k];
        }

        err = getforce_tracers(&xt3[0], &vt3[0], acc3, gal, ngals);
        for (int k=0; k<3; k++) {                /* Runge-Kutta formula */
            (*test_particles).pos[3*n+k] += dt/6.0*(v[k]+2.0*(vt1[k]+vt2[k])+vt3[k]);
            (*test_particles).vel[3*n+k] += dt/6.0*(acc0[k]+2.0*(acc1[k]+acc2[k])+acc3[k]);
        }
    }
    return err;
}

int getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, int ngals){

    int i;
    int err = 0;
    double r, rx, ry, rz;
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double vr = 0.0;
    // get the force from all other galaxies
    for (i=0; i<ngals; i++){
        if (i != gal_num){  // skip itself
            //Hernquist bulge
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
                 pow(*(x+1) - gal[i].pos[1], 2) +
                 pow(*(x+2) - gal[i].pos[2], 2));
            if (r > RMIN){ // control for very close encounters
                rx = (*x - gal[i].pos[0]);
                ry = (*(x+1) - gal[i].pos[1]);
                rz = (*(x+2) - gal[i].pos[2]);
            } else {
                printf("Very close\n");
                rx = r;
                ry = r;
                rz = r;
            }

    	    ax += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
    	        (r+gal[i].b1_LMJ))*rx/r;
    	    ay += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
    	        (r+gal[i].b1_LMJ))*ry/r;
    	    az += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
    	        (r+gal[i].b1_LMJ))*rz/r;

            //Miyamato disk
    	    r = sqrt(pow(*x - gal[i].pos[0], 2) +
    	             pow(*(x+1) - gal[i].pos[1], 2) +
    	             pow(gal[i].a2_LMJ + sqrt(pow(*(x+2) - gal[i].pos[2], 2)
    	                                      + pow(gal[i].b2_LMJ, 2)), 2)
    	             );
            if (r > RMIN){ // control for very close encounters

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
            if (r > RMIN){ // control for very close encounters
                rx = (*x - gal[i].pos[0]);
                ry = (*(x+1) - gal[i].pos[1]);
                rz = (*(x+2) - gal[i].pos[2]);
            } else {
                rx = r;
                ry = r;
                rz = r;
            }
            ax += halo_acc(r, gal[i], rx, 0.0);  //x vector
            ay += halo_acc(r, gal[i], ry, 0.0);  // y vector
            az += halo_acc(r, gal[i], rz, 0.0);  // z vector
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
                if (err) {
                    return err;
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

int getforce_tracers(double *x, double *v, double *a, struct Gal *gal, int ngals){

    int i;
    int err = 0;
    double r, rx, ry, rz;
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
    // get the force from all other galaxies
    for (i=0; i<ngals; i++){
        //Hernquist bulge
        r = sqrt(pow(*x - gal[i].pos[0], 2) +
             pow(*(x+1) - gal[i].pos[1], 2) +
             pow(*(x+2) - gal[i].pos[2], 2));
        if (r > RMIN){ // control for very close encounters
            rx = (*x - gal[i].pos[0]);
            ry = (*(x+1) - gal[i].pos[1]);
            rz = (*(x+2) - gal[i].pos[2]);
        } else {
            rx = r;
            ry = r;
            rz = r;
        }

        ax += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
            (r+gal[i].b1_LMJ))*rx/r;
        ay += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
            (r+gal[i].b1_LMJ))*ry/r;
        az += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
            (r+gal[i].b1_LMJ))*rz/r;

        //Miyamato disk
        r = sqrt(pow(*x - gal[i].pos[0], 2) +
                 pow(*(x+1) - gal[i].pos[1], 2) +
                 pow(gal[i].a2_LMJ + sqrt(pow(*(x+2) - gal[i].pos[2], 2)
                                          + pow(gal[i].b2_LMJ, 2)), 2)
                 );
        if (r > RMIN){ // control for very close encounters

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
        if (r > RMIN){ // control for very close encounters
            rx = (*x - gal[i].pos[0]);
            ry = (*(x+1) - gal[i].pos[1]);
            rz = (*(x+2) - gal[i].pos[2]);
        } else {
            rx = r;
            ry = r;
            rz = r;
        }
        ax += halo_acc(r, gal[i], rx, 0.0);  //x vector
        ay += halo_acc(r, gal[i], ry, 0.0);  // y vector
        az += halo_acc(r, gal[i], rz, 0.0);  // z vector
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
    double sigma,density;
    double dyn_L, dyn_C, dyn_alpha;
    double mhalo = gal.mhalo;

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
            dyn_L = gal.dyn_L_eq;
            dyn_C = gal.dyn_C_eq;
            dyn_alpha = gal.dyn_alpha_eq;
        } else {
            dyn_L = gal.dyn_L_uneq;
            dyn_C = gal.dyn_C_uneq;
            dyn_alpha = gal.dyn_alpha_uneq;
        }
        coulomb = fmax(dyn_L, pow(log(r/(dyn_C*r_gal)),
                                    dyn_alpha));
    }
    /*
    * X and Velocity dispersion
    */

    /*if (gal.tidal_trunc){
        sigma = halo_sigma_trunc(r, gal);
    } else {
        sigma = halo_sigma(r, gal);
    }
    */
    sigma = halo_sigma(r, gal);

    if (sigma < 0){
        printf("sigma: %g\n", sigma);
        return 1;
    }

    sigma = 3.0*sqrt(sigma);

    //halo_sigma_old(r, gal, &sigma);
    double X = vr/(sqrt(2.0)*sigma);

    density = halo_density(r, gal);

    *ax += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vx/pow(vr, 3);

    *ay += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vy/pow(vr, 3);

    *az += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vz/pow(vr, 3);

    return 0;
}


double tidal_condition(double x, void *params)
{
    double *p = (double *)params;
    double MD, MG;
  // p[0] = MD, p[1] = MG, p[2] = aD, p[3] = aG, p[4] = r, p[5] = gammaD, p[6] = gammaG
    MD = p[0]*pow(x/(x+p[2]), 3-p[5]);

    MG = p[1]*pow((p[4]-x)/(p[4]-x+p[3]), 3-p[6]);

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
    double p[7] = {galD.mhalo, galG.mhalo, galD.r_halo, galG.r_halo, r,
                   galD.gamma, galG.gamma};
    F.function = &tidal_condition;
    F.params = (void *)p;
    gsl_set_error_handler(&custom_gsl_error_handler);

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    status = gsl_min_fminimizer_set(s, &F, m, a, b);
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
      return -2.0; // GSL reports non-zero error codes
    }
    gsl_min_fminimizer_free (s);

    return m;
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
        getforce_gals(gal[n].pos, gal[n].vel, acc0, n, gal, ngals);
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


void record_snapshot(int ngals, struct Gal *gal, double t, int snapnumber, struct Snapshot **output_snapshot){
    int n, i;
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

void write_tracers(struct Params parameters, struct Gal *gal, double t, int snapnumber){
    int ngals = parameters.ngals;
    char *folder = parameters.outputdir;
    FILE *snapfile;
    char snapname[50];
    sprintf(snapname, "%ssnapshot.tracers.csv.%03d", folder, snapnumber);
    snapfile = fopen(snapname, "w");
    fprintf(snapfile,"NAME,X,Y,Z,VX,VY,VZ,T\n");
    for (int n=0; n<ngals; n++){
        if (gal[n].test_particles.nparticles > 0){
            for (int m=0; m<gal[n].test_particles.nparticles; m++){
                fprintf(snapfile,"%s,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f\n",
                        gal[n].name,
                        gal[n].test_particles.pos[m*3+0],
                        gal[n].test_particles.pos[m*3+1],
                        gal[n].test_particles.pos[m*3+2],
                        gal[n].test_particles.vel[m*3+0],
                        gal[n].test_particles.vel[m*3+1],
                        gal[n].test_particles.vel[m*3+2],
                        t);
            }
        }
    }
    fclose(snapfile);
}
