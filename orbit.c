#include "orbit.h"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int orbit(int int_mode,
          int ngals,
          struct Params parameters,
          struct Gal *gal,
          double* output_pos,
          double* output_vel){


    int ratio, n, i;
    double tpast = parameters.tpast;  // add to orbit.h
    double tfuture = parameters.tfuture; // add to orbit.h
    double dt0 = parameters.dt0; // add to orbit.h
    double dtout = parameters.dtout;
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;

    //get position of cluster at t = -tpast
    double sign, tmax, dtoutt, t;
    int err; 
    if (tpast < 0.0) { 
        sign = -1.0;
        tmax = tpast;
        // If we don't want to integrate forward, then save snapshots as we go back
        if (tfuture <= 0.0){ 
            dtoutt = -1.0*dtout;
        } else {  // otherwise just save a snapshot at the end of the backward integration
            dtoutt = tpast;
        }
        t = tstart;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal, parameters, sign);
        printf("%d\n", err);
    }
    if (tfuture > 0.0) {
    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
        sign = 1.0;
        dtoutt = dtout;
        tmax = tfuture;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal, parameters, sign);
        printf("%d\n", err);
    }
    /*}
    else
    {
        double sign;
        if (tpast < 0.0)
            sign = -1.0;
        else
            sign = 1.0;
        double tmax = tpast;
        double dtoutt = dtout;
        double t = tstart;
        int err;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal, parameters, sign);
        printf("%d", err);
    }*/

       //}
    for (n=0; n<ngals; n++){
        for (i=0; i<3; i++){
            output_pos[i+n*3] = gal[n].pos[i];
            output_vel[i+n*3] = gal[n].vel[i];
        }
    }
    free(gal);
    return 0;
}


/* --------------- extrapolation method --------------- */
int rk4_drv(double *t,
            double tmax,
            double dtout,
            double dt0,
            double mdiff,
            struct Gal *gal,
            struct Params parameters,
            double sign){
        struct OrbitStats stats;
        int snapnum = 0;
	double tout, diff, dt = 0.0;
	double xe1[3], ve1[3], difftemp, dist;
        double old_dist = -1.0;
        //double rt = 1e+5;
        double rt_temp, r;
	int dir = -1;
        int k, n, m;
	int ngals = parameters.ngals;
        //char name[50];
        //sprintf(name, "temp.dat");
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
		    for (k=0;k<3;k++) {
                        (*(gal+n)).post[k] = (*(gal+n)).pos[k];
                        (*(gal+n)).velt[k] = (*(gal+n)).vel[k];
			xe1[k]=(*(gal+n)).pos[k];
			ve1[k]=(*(gal+n)).vel[k];
		    }
                    if (VARIABLE_TIMESTEPS) {
		        do_step(dt, xe1, ve1, n, gal, parameters);      /* One full step */
		        do_step(0.5*dt, (*(gal+n)).post, (*(gal+n)).velt, n, gal, parameters);  /* Two half steps */
		        do_step(0.5*dt, (*(gal+n)).post, (*(gal+n)).velt, n, gal, parameters);
		        difftemp = sqrt(pow(xe1[0] - (*(gal+n)).post[0],2) +
		                        pow(xe1[1] - (*(gal+n)).post[1],2) +
		                        pow(xe1[2] - (*(gal+n)).post[2],2));
                        if (difftemp > diff) {diff = difftemp;} // hold highest value to compare with mdiff below
                    } else {
		        do_step(dt, (*(gal+n)).post, (*(gal+n)).velt, n, gal, parameters);      /* One full step */
                    }
                }
                // end loop over each galaxy
		if (!VARIABLE_TIMESTEPS || (diff<=mdiff)) {         /* Is difference below accuracy threshold? */
		    *t+=dt;/* If yes -> continue and double step size */
		    //TODO: update the test particles here
                    for (n=0; n<ngals; n++){
		        for (k=0;k<3;k++) {
			    (*(gal+n)).pos[k]=(*(gal+n)).post[k];
			    (*(gal+n)).vel[k]=(*(gal+n)).velt[k];
		        }
		    }
                    // Tidal stripping
                    //if (dt > 0) {
                    for  (n=0; n<ngals; n++){
                        if (gal[n].tidal_trunc == 1) {// tidal truncation turned on
                            for (m=0; m<ngals; m++){ // look for all galaxies with dynamic friction turned on
                                if ((m != n) && (gal[m].dyn_fric == 1)){ 
                                    r = sqrt(pow(gal[n].pos[0]-gal[m].pos[0], 2) +
                                             pow(gal[n].pos[1]-gal[m].pos[1], 2) +
                                             pow(gal[n].pos[2]-gal[m].pos[2], 2));
                                    rt_temp = calc_rt(r, fmin(gal[n].rt, gal[n].r_halo), gal[m], gal[n]);
                                    if (sign > 0){  // integrating forward
                                        gal[n].rt = fmin(gal[n].rt, rt_temp);
                                    } else if (sign < 0){  // integrating backward
                                        gal[n].rt = fmax(gal[n].rt, rt_temp);
                                    }
                                    if (gal[n].halo_type == 1){ // Dehnen
                                        gal[n].mhalo = gal[n].minit*pow(gal[n].rt/(gal[n].rt+gal[n].r_halo), 3-gal[n].gamma); // Dehnen
                                    } else if (gal[n].halo_type == 2){ // NFW
                                        gal[n].mhalo = gal[n].minit*(log(1+(gal[n].rt)/gal[n].r_halo)-gal[n].rt/(gal[n].r_halo+gal[n].rt))/(log(1+gal[n].c_halo) -gal[n].c_halo/(1+gal[n].c_halo));
                                    } else { // Plummer
                                        gal[n].mhalo = gal[n].minit*pow(gal[n].rt, 3)/pow(pow(gal[n].r_halo, 2) + pow(gal[n].rt, 2), 1.5);
                                    }
                                }
                            }
                        }
                    }
                    //}
		    if (VARIABLE_TIMESTEPS) { //If we are with the threshold double the timestep and continue
                        dt = dt*2.0;
                    }
		} else {
		    dt = dt/2.0;  // if we are outside the threshold halve the timestep and try again
		}
                // Abort if the timestep ever gets too low
		if (sign*dt < 0.01*dt0 && !laststep) {
		    printf("Aborted... dt = %lf (>%lf, %lf)\n", dt, dt0, sign);
		    return 1;
		}
		count++;
                // round to the end of simulation time
		if (sign*dt > dtmax) {
		    dt = sign*dtmax;
		}

		} while (diff>mdiff);       /* Go through loop once and only repeat if difference is too large */
	    if (sign**t>=sign*(tout)) {
                for (n=0; n<ngals; n++) {
                    if (gal[n].tidal_trunc == 1) {
                        printf("rt: %10.5f, m: %10.5f, t: %10.5f\n", gal[n].rt, gal[n].mhalo, *t);
                    }
                }
                if (snapshot){
                    write_snapshot(parameters, gal, *t, snapnum);
                    snapnum += 1;
                }
                if (orbit_stats) {
                    //record orbit stats here (pericenters and apocenters etc.)
                    dist = sqrt( pow((*(gal+ref_gal)).pos[0]-(*(gal+test_gal)).pos[0], 2) +
                                 pow((*(gal+ref_gal)).pos[1]-(*(gal+test_gal)).pos[1], 2) +
                                 pow((*(gal+ref_gal)).pos[2]-(*(gal+test_gal)).pos[2], 2));
                    if (old_dist >= 0){
                        if (dir >= 0){
                        // already have a direction set
                            if ((dir = 1) && (dist >= old_dist)) {
                            // moving in and got further away
                                dir = 0;
                                stats.dir = dir;
                                stats.pericenters += 1;
                            } else if ((dir = 0) && (dist <= old_dist)) {
                            // moving out and got closer
                                dir = 1;
                                stats.apocenters += 1;
                                stats.dir = dir;
                            }
                        } else {
                        // set direction
                            if (dist < old_dist) {
                                dir = 1;  // inward
                                stats.dir = dir;
                            } else {
                                dir = 0;  // outward
                                stats.dir = dir;
                            }
                        }
                    }
                    old_dist = dist;
                }
                tout+=dtout;              /* increase time of output/next insertion */
            }


        } while (sign**t<sign*(tmax));
        // write final snapshot
        if (snapshot) write_snapshot(parameters, gal, *t, snapnum);
	return 0;

}


/* ---------- advancement ---------- */
void do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, struct Params parameters) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;

    hh = dt*0.5;
    if (RK4) {
    	getforce_gals(x, v, acc0, gal_num, gal, parameters);
    	for (k=0;k<3;k++) {                /* first half-step */
    	    xt1[k] = *(x+k)+hh**(v+k);
    	    vt1[k] = *(v+k)+hh**(acc0+k);
    	}

    	getforce_gals(&xt1[0], &vt1[0], acc1, gal_num, gal, parameters);
    	for (k=0;k<3;k++) {                /* second half-step */
    	    xt2[k] = *(x+k)+hh*vt1[k];
    	    vt2[k] = *(v+k)+hh**(acc1+k);
    	}

    	getforce_gals(&xt2[0], &vt2[0], acc2, gal_num, gal, parameters);
    	for (k=0;k<3;k++) {                /* third half-step with results of second half-step */
    	    xt3[k] = *(x+k)+dt*vt2[k];
    	    vt3[k] = *(v+k)+dt**(acc2+k);
    	}

	getforce_gals(&xt3[0], &vt3[0], acc3, gal_num, gal, parameters);
	for (k=0;k<3;k++) {                /* Runge-Kutta formula */
	    *(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
	    *(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
	}
    }
    else {  // modified leapfrog
        // ai
        getforce_gals(x, v, acc0, gal_num, gal, parameters);
        // vi+1/2 and xi+1/2
        for (k=0;k<3;k++) {
            vt1[k] = *(v+k)+ *(acc0+k)*hh;
            xt1[k] = *(x+k)+ *(v+k)*hh;
        }
        // ai+1/2
        getforce_gals(xt1, vt1, acc1, gal_num, gal, parameters);
        //vi+1 and xi+1
        for (k=0;k<3;k++) {
            vt1[k] = *(v+k)+ *(acc1+k)*dt;
            *(x+k) += 0.5*(*(v+k)+vt1[k])*dt;
            *(v+k) = vt1[k];
        }
    }
}


void getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, struct Params parameters){
    getforce(x, v, a, parameters, gal[gal_num]);
    int i;
    double r;
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double vr = 0.0;
    double constant;

    int ngals = parameters.ngals;
    // get the force from all other galaxies
    for (i=0; i<ngals; i++){
        if (i != gal_num){  // skip itself
            //Hernquist bulge
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
                 pow(*(x+1) - gal[i].pos[1], 2) +
                 pow(*(x+2) - gal[i].pos[2], 2));

	    ax += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
	        (r+gal[i].b1_LMJ))*(*x - gal[i].pos[0])/r;
	    ay += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
	        (r+gal[i].b1_LMJ))*(*(x+1) - gal[i].pos[1])/r;
	    az += -G*gal[i].M1_LMJ/((r+gal[i].b1_LMJ)*
	        (r+gal[i].b1_LMJ))*(*(x+2) - gal[i].pos[2])/r;
          //  if (r == r) {
          //      printf("Bulge %10.5f %10.5f %10.5f %10.5f\n", r, ax, ay, az); 
          //  }
            //Miyamato disk
	    r = sqrt(pow(*x - gal[i].pos[0], 2) +
	             pow(*(x+1) - gal[i].pos[1], 2) +
	             pow(gal[i].a2_LMJ + sqrt(pow(*(x+2) - gal[i].pos[2], 2)
	                                      + pow(gal[i].b2_LMJ, 2)), 2)
	             );

	    ax += -G*gal[i].M2_LMJ/(r*r*r) * (*x - gal[i].pos[0]);
	    ay += -G*gal[i].M2_LMJ/(r*r*r) * (*(x+1) - gal[i].pos[1]);
	    az += -G*gal[i].M2_LMJ/(r*r*r) * (gal[i].a2_LMJ + sqrt(pow(*(x+2) - gal[i].pos[2], 2)
	                                                    + pow(gal[i].b2_LMJ, 2)))
	                                             / sqrt(pow(*(x+2) - gal[i].pos[2], 2)
	                                                    + pow(gal[i].b2_LMJ, 2))
	                                             * (*(x+2) - gal[i].pos[2]);
          //  if (r == r) {
          //      printf("disk %10.5f %10.5f %10.5f %10.5f\n", r, ax, ay, az); 
          //  }
            // Dark Matter Halo
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
               pow(*(x+1) - gal[i].pos[1], 2) +
               pow(*(x+2) - gal[i].pos[2], 2));

            if (gal[i].halo_type == 1) { // Dehnen
                ax += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*x - gal[i].pos[0]);
                ay += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*(x+1) - gal[i].pos[1]);
                az += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*(x+2) - gal[i].pos[2]);
            } else if (gal[i].halo_type == 2) {  // NFW
                constant = -G*gal[i].mhalo/
                        (log(1+gal[i].c_halo)-gal[i].c_halo/(1+gal[i].c_halo));
                ax += constant * (log(1.0 + r/gal[i].r_halo)/r -
                                1.0/(gal[i].r_halo+r)) * (*x - gal[i].pos[0])/pow(r, 2);
                ay += constant * (log(1.0 + r/gal[i].r_halo)/r -
                                1.0/(gal[i].r_halo+r)) * (*(x+1) - gal[i].pos[1])/pow(r, 2);
                az += constant * (log(1.0 + r/gal[i].r_halo)/r -
                                1.0/(gal[i].r_halo+r)) * (*(x+2) - gal[i].pos[2])/pow(r, 2);

            } else { // Plummer sphere
                ax += -2.0*G*gal[i].mhalo* (*x - gal[i].pos[0])/
                    pow(pow(gal[i].r_halo, 2)+pow(r, 2), 1.5);
                ay += -2.0*G*gal[i].mhalo* (*(x+1) - gal[i].pos[1])/
                    pow(pow(gal[i].r_halo, 2)+pow(r, 2), 1.5);
                az += -2.0*G*gal[i].mhalo* (*(x+2) - gal[i].pos[2])/
                    pow(pow(gal[i].r_halo, 2)+pow(r, 2), 1.5);
            }
          //  if (r == r) {
            //    printf("Halo %10.5f %10.5f %10.5f %10.5f\n", r, ax, ay, az); 
           // }
            // dynamical friction
            if (gal[i].dyn_fric == 1) {// is dynamical friction turned on for this galaxy?
                //relative velocity
                vx = (*v - gal[i].vel[0]);
                vy = (*(v+1) - gal[i].vel[1]);
                vz = (*(v+2) - gal[i].vel[2]);
                vr = sqrt(vx*vx + vy*vy + vz*vz);

                dynamical_friction(r, vx, vy, vz, vr,
                                   &ax, &ay, &az,
                                   gal[i].halo_type, gal[i].mhalo,
                                   gal[i].r_halo, gal[i].gamma, gal[i].c_halo,
                                   gal[gal_num].mhalo, gal[gal_num].r_halo);
            
                   
            }
        }
    }
    // update acceleration
    *(a+0) += ax;
    *(a+1) += ay;
    *(a+2) += az;

}


/* ----------- force ----------- */
void getforce(double *x, double *v, double *a, struct Params parameters, struct Gal gal){
	double r1, r2, r3;
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;

    // Potential for MW
	//Hernquist bulge
    r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));

    ax += -G*parameters.M1_LMJ/((r1+parameters.b1_LMJ)*
            (r1+parameters.b1_LMJ))**x/r1;
    ay += -G*parameters.M1_LMJ/((r1+parameters.b1_LMJ)*
            (r1+parameters.b1_LMJ))**(x+1)/r1;
    az += -G*parameters.M1_LMJ/((r1+parameters.b1_LMJ)*
            (r1+parameters.b1_LMJ))**(x+2)/r1;

    //Miyamato disk
    r2 = sqrt(pow(*x, 2) +
              pow(*(x+1), 2) +
              pow(parameters.a2_LMJ + sqrt(pow(*(x+2), 2) + pow(parameters.b2_LMJ, 2)), 2));

    ax += -G*parameters.M2_LMJ/(r2*r2*r2) * *x;
    ay += -G*parameters.M2_LMJ/(r2*r2*r2) * *(x+1);
    az += -G*parameters.M2_LMJ/(r2*r2*r2) * (parameters.a2_LMJ + sqrt(pow(*(x+2), 2)
                + pow(parameters.b2_LMJ, 2)))/
            sqrt(pow(*(x+2), 2) + pow(parameters.b2_LMJ, 2)) * *(x+2);

    // Dark matter halo
    if (parameters.Mhalo > 0.0){ // Avoid divide by zero errors

        r3 = sqrt(pow(*x, 2) + pow(*(x+1), 2) + pow(*(x+2)/parameters.q_halo, 2));

        if (parameters.halo_type == 2) {  //NFW
            double constant = -G*parameters.Mhalo/
                        (log(1+parameters.c_halo)-parameters.c_halo/(1+parameters.c_halo));
            ax += constant * (log(1.0 + r3/parameters.r_halo)/r3 -
                                1.0/(parameters.r_halo+r3)) * *x/pow(r3, 2);
            ay += constant * (log(1.0 + r3/parameters.r_halo)/r3 -
                                1.0/(parameters.r_halo+r3)) * *(x+1)/pow(r3, 2);
            az += constant * (log(1.0 + r3/parameters.r_halo)/r3 -
                                1.0/(parameters.r_halo+r3)) * *(x+2)/
                                pow(parameters.q_halo*parameters.q_halo*r3, 2);

        } else { // Dehnen
            ax += -G*parameters.Mhalo * pow(r3/(parameters.r_halo + r3), -parameters.gamma)/
                                    pow(parameters.r_halo + r3, 3) * *x;
            ay += -G*parameters.Mhalo * pow(r3/(parameters.r_halo + r3), -parameters.gamma)/
                                    pow(parameters.r_halo + r3, 3) * *(x+1);
            az += -G*parameters.Mhalo * pow(r3/(parameters.r_halo + r3), -parameters.gamma)/
                                    pow(parameters.r_halo + r3, 3) * *(x+2)/pow(parameters.q_halo, 2);
        }
        if (DYNAMICALFRICTION_MAIN) {
            //relative velocity
            double vr = sqrt(*v* *v + *(v+1)* *(v+1) + *(v+2) * *(v+2));
            dynamical_friction(r3, *v, *(v+1), *(v+2), vr,
                               &ax, &ay, &az,
                               parameters.halo_type, parameters.Mhalo,
                               parameters.r_halo, parameters.gamma, parameters.c_halo,
                               gal.mhalo, gal.r_halo);


        }
    }
    *(a+0) = ax;
    *(a+1) = ay;
    *(a+2) = az;
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


// Calculate Dynamical Friction acceleration
void dynamical_friction(double r, double vx, double vy, double vz, double vr,  // orbit velocity and radius
                        double *ax, double *ay, double *az,  // accelerations update in function
                        int halo_type, double mhalo, double r_halo, double gamma, double c_halo, // Halo properties
                        double m_gal, double r_gal){  // companion mass and scale length
    double sigma = 0.0;
    double density = 0.0;
    double dyn_L, dyn_C, dyn_alpha;
    /*
    * XXXXXXXXXXXXXXXXX
    * COULOMB LOGARITHM
    * XXXXXXXXXXXXXXXXX
    */
    double coulomb;
    // alternative methods of coulomb logarithm
    if (mhalo/m_gal < 0.2){
        coulomb = 0.0;  // don't have a small galaxy act on a big one
    } else {
        // van den Marel et al. 2012 eq. A1 and discussion in Appendix A
        if (abs(mhalo/m_gal - 1.0) < 0.3){ // within 30% of each other
            //dyn_L = 0.02;
            dyn_L = 0.0;
            dyn_C = 1.0;
            //dyn_C = 0.17;
           // dyn_alpha = 0.15;
            dyn_alpha = 2.5;
        } else {
            //Test against gadget
            //dyn_L = 0.0;
            //dyn_C = 0.1;
            //dyn_alpha = 3.5;
            //van der Marel 2012
            //dyn_C =  1.22;
            //dyn_alpha = 1.0;
            //K+13 below
            dyn_L = 0.0;
            dyn_C = 1.6*3.0/r_gal;  // based off K+13 etc.
            dyn_alpha = 1.0;
        }
        coulomb = fmax(dyn_L, pow(log(r/(dyn_C*r_gal)),
                                    dyn_alpha));
    }
    /*
    * XXXXXXXXXXXXXXXXXXXXXXXXX
    * X and Velocity dispersion
    * XXXXXXXXXXXXXXXXXXXXXXXXX
    */
    if (halo_type == 1){
    // Hernquist 1D velocity dispersion - CForm from Mathematica -- only good for gamma == 1
        sigma = 3.0*sqrt(G*mhalo*r*pow(r_halo + r, 3)*
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
        sigma = 3.0* VMAX *1.4393*pow(r, 0.354)/(1+1.1756*pow(r, 0.725));
    } else { //Plummer
        sigma = 3.0*sqrt(pow(r_halo, 5)*G*mhalo*
                    pow(1+pow(r/r_halo, 2), 2.5)/
                    (6.0*pow(pow(r_halo, 2)+pow(r, 2), 3)));
    }

    double X = vr/(sqrt(2.0)*sigma);
    if (halo_type == 1) { //Hernquist
        density = (3.0 - gamma)*mhalo/(4.0*Pi)*
              r_halo/(pow(r, gamma)*pow(r + r_halo, 4-gamma));
    } else if (halo_type == 2) { //NFW
        // where rho0 = Mvir/(4*Pi*a^3)/(log(1+c)-c/(1+c))
        density = mhalo/(4.0*Pi*pow(r_halo, 3))/
                        (log(1+c_halo)-c_halo/(1+c_halo))/
                        (r/r_halo*pow(1+r/r_halo, 2));
    } else { //Plummer
        density = 3*mhalo/(4*Pi*pow(r_halo, 3)*
                pow(1+pow(r/r_halo, 2), 2.5));
    }
    *ax += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vx/pow(vr, 3);

    *ay += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vy/pow(vr, 3);

    *az += -4.0*Pi*G*G*m_gal*density*coulomb*
        (erf(X) - 2.0*X/sqrt(Pi)*exp(-X*X))*vz/pow(vr, 3);

}

double tidal_condition (double x, void * params)
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

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status 
        = gsl_min_test_interval (a, b, 0.001, 0.001);

    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);

  return m;
}
