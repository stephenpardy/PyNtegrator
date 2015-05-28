#include "orbit.h"
#include "orbit_utils.h"
//#include <Python.h>
//#include <numpy/arrayobject.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int orbit(double mass_cluster,
          double pm_mu_delta,
          double pm_mu_alphacosdelta,
          double distance_cluster,
          double q_halo,
          double r_halo,
          double tpast,
          double sigma_x,
          double sigma_v,
          double sigma_vx,
          double sigma_mu) {

    //round simulation time to multiple of output time
    struct Gal gal;
    int ratio;
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;

    //put the MCMC parameters in the appropriate array for passing it to the integrator
    double parameter[17];
    parameter[0] = mass_cluster;            //mass of cluster
    parameter[1] = pm_mu_delta;             //delta component of proper motion
    parameter[2] = pm_mu_alphacosdelta;     //alpha component of proper motion
    parameter[3] = distance_cluster;        //distance cluster-sun
    parameter[4] = q_halo;                 //flattening of dark halo
    parameter[5] = r_halo;                 //concentration of dark halo
    parameter[6] = tpast;                //integration time
    parameter[7] = sigma_x;                //smoothing position of overdensities
    parameter[8] = sigma_v;                //smoothing v_r
    parameter[9] = sigma_vx;               //smoothing position of v_r stars
    parameter[10] = sigma_mu;              //smoothing proper motions of v_r stars


    //convert coordinates into Cartesian frame
    double dsuntemp;
    double vrsuntemp;
    double vrtemp;
    double ltemp;
    double btemp;
    double lcosbtemp;
    double RAtemp;
    double DECtemp;
    double mu_alphatemp;
    double mu_alphacosdeltatemp;
    double mu_deltatemp;
    double mutemp;
    double PAtemp;
    double vLSRtemp;
	dsuntemp = distance_cluster*1000.;
	mu_alphatemp = 0.0;
	mu_alphacosdeltatemp = pm_mu_alphacosdelta*1000.0;
	mu_deltatemp = pm_mu_delta*1000.0;
	vrsuntemp = vLSR;
	ltemp = l;
	btemp = b;
	lcosbtemp = lcosb;
	vLSRtemp = vLSR;
	convert(gal.pos,
            gal.vel,
            &dsuntemp,
            &vrsuntemp,
            &vrtemp,
            &ltemp,
            &btemp,
            &lcosbtemp,
            &RAtemp,
            &DECtemp,
            &mu_alphatemp,
            &mu_alphacosdeltatemp,
            &mu_deltatemp,
            &mutemp,
            &PAtemp,
            2,
            2,
            0,
            vLSRtemp,
            rgalsun);


    //get the acceleration at cluster's position
	double atemp[3];
	getforce(gal.pos, gal.vel, atemp, parameter);


    //get position of cluster at t = -tpast
    double sign = -1.0;
    double tmax = tpast;
    double dtoutt = tpast;
    double t = tstart;
    rk4_drv(&t, tmax, dtoutt, mdiff, &gal, 1, sign, parameter);

    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
    sign = 1.0;
    dtoutt = dtout;
    tmax = tfuture;
    int err;
    err = rk4_drv(&t, tmax, dtoutt, mdiff, &gal, 1, sign, parameter);

    //}

	return 0;
}

/* --------------- extrapolation method --------------- */
int rk4_drv(double *t,
            double tmax,
            double dtout,
            double mdiff,
            struct Gal *gal,
            int ngals,
            double sign,
            double *parameter){

	double tout, diff, dt = 0.0;
	double atemp[3],difftemp;
	xe1[3], xe2[3], ve1[3], ve2[3]
	int k, n;
        char name[50];
        sprintf(name, "temp.dat");
	//initialize output/insertion of tail particles
	tout = *t;		    /* time of next output/insertion */
	dt = sign*dt0;                /* initial time step */
	//integrate cluster
	do {

		/***********
		 * CLUSTER *
		 ***********/
            //advance cluster particle
	    int count = 0;
            int laststep = 0;
	    do {
	        //if (sign*(*t+dt) > sign*tout) {
		//	dt = tout-*t;
                //    laststep = 1;
                //    printf("setting dt: %lf \n", dt);
			//}
		difftemp = 0.0;
		diff = 0.0;
                // loop over each galaxy
                for (n=0; n<ngals; n++){
		    for (k=0;k<3;k++) {
			xe1[n][k]=(*(gal+n)).pos[k];
			xe2[n][k]=(*(gal+n)).pos[k];
			ve1[n][k]=(*(gal+n)).vel[k];
			ve2[n][k]=(*(gal+n)).vel[k];
		    }
		    do_step(dt,xe1[n],ve1[n],parameter);      /* One full step */

		    do_step(0.5*dt, xe2[n], ve2[n], parameter);  /* Two half steps */
		    do_step(0.5*dt, xe2[n], ve2[n], parameter);
		    difftemp = sqrt(pow(xe1[n][0] - xe2[n][0],2) + pow(xe1[n][1] - xe2[n][1],2) + pow(xe1[n][2] - xe2[n][2],2));
                    if (difftemp > diff) {diff = difftemp;} // hold highest value to compare with mdiff below
                }
                // end loop over each galaxy
		if (diff<mdiff) {         /* Is difference below accuracy threshold? */
		    *t+=dt;
                    for (n=0; n<ngals; n++){
		        for (k=0;k<3;k++) {
			    (*(gal+n)).pos[k]=xe2[n][k];          /* If yes -> continue and double step size */
			    (*(gal+n)).vel[k]=ve2[n][k];
		        }
		    }
		    dt = dt*2.0;
		} else {
		    dt = dt/2.0;
		}

		if (sign*dt < 0.01*dt0 && !laststep) {
		    printf("Aborted... dt = %lf (>%lf, %lf)\n", dt, dt0, sign);
		    return 1;
		}
		count++;

		if (sign*dt > dtmax) {
		    dt = sign*dtmax;
		}

		} while (diff>mdiff);       /* Go through loop once and only repeat if difference is too large */
            if ((sqrt(pow((*(gal+n)).pos[0],2)+pow((*(gal+n)).pos[1],2)+pow((*(gal+n)).pos[2],2)) < Rgalmin) ||
                (sqrt(pow((*(gal+n)).pos[0],2)+pow((*(gal+n)).pos[1],2)+pow((*(gal+n)).pos[2],2)) > Rgalmax)  ) return 1;  /* Abort if r is too extreme */
	    if (sign**t>=sign*(tout)) {
                // save info on both gals
                printf("x,y,z, t = %lf, %lf, %lf %lf \n", (*(gal+n)).pos[0], (*(gal+n)).pos[1], (*(gal+n)).pos[2], *t);
                tout+=dtout;              /* increase time of output/next insertion */
            }

        } while (sign**t<sign*(tmax));


	//last printout
	//double dsuntemp, vrsuntemp, vrtemp, ltemp, btemp, lcosbtemp, RAtemp, DECtemp, mu_alphatemp, mu_alphacosdeltatemp, mu_deltatemp, mutemp, PAtemp, vLSRtemp;
	//vLSRtemp = parameter[7];
	//getforce((*gal).pos,(*gal).vel,atemp,parameter);
	//convert((*gal).pos, (*gal).vel, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp,rgalsun);
	//printf("\n%.2f\n\n", sqrt(atemp[0]*atemp[0]+atemp[1]*atemp[1]+atemp[2]*atemp[2]));

	return 0;

}


/* ----------- force ----------- */
void getforce(double *x, double *v, double *a, double *parameter){
	double r1, r2, r3;
	double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;

    // Potential for MW
	//Hernquist bulge
	r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));

	a1x = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**x/r1;
	a1y = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+1)/r1;
	a1z = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+2)/r1;

	//Miyamato disk
	r2 = sqrt(*x * *x + *(x+1) * *(x+1) + (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))*(a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ)));

	a2x = -G*M2_LMJ/(r2*r2*r2) * *x;
	a2y = -G*M2_LMJ/(r2*r2*r2) * *(x+1);
	a2z = -G*M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))/sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ) * *(x+2);

	//NFW Halo
	r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[4] * *(x+2)/parameter[4]);

    a3x = -G*Mhalo/r3 * (log(1.0 + r3/parameter[5])/r3 - 1.0/(parameter[5]+r3)) * *x/r3;
    a3y = -G*Mhalo/r3 * (log(1.0 + r3/parameter[5])/r3 - 1.0/(parameter[5]+r3)) * *(x+1)/r3;
    a3z = -G*Mhalo/r3 * (log(1.0 + r3/parameter[5])/r3 - 1.0/(parameter[5]+r3)) * *(x+2)/(parameter[4]*parameter[4]*r3);

	*(a+0) = a1x + a2x + a3x;
	*(a+1) = a1y + a2y + a3y;
	*(a+2) = a1z + a2z + a3z;

}

/*
void getforce_gals(double *x, double *v, double *a, double *parameter){
    double r1, r2, r3;
    double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
    double q1, q2, qz, phi, C1, C2, C3, rhalo, vhalo;
//    double r200, c, deltachar;
//    double k = 1.3e-7;
    //NFW Halo
    r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);

    a3x = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *x/r3;
    a3y = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+1)/r3;
    a3z = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+2)/(parameter[10]*parameter[10]*r3);

    *(a+0) = a3x;
    *(a+1) = a3y;
    *(a+2) = a3z;

}
*/

/* ----------- force tail -----------
void getforce_testparticle(double *x, double *v, double *a, double *xc, double *parameter){
    double r4;
    double a4x, a4y, a4z;
	double actualmass;
	actualmass = parameter[0]+parameter[9]*parameter[8];
    // MW
    getforce(x, v, a, parameter);

	//Galaxies
    //getforce_gals(x, v, a4, parameter);
	r4 = sqrt((*x-*xc) * (*x-*xc) +
              (*(x+1)-*(xc+1)) * (*(x+1)-*(xc+1)) +
              (*(x+2)-*(xc+2)) * (*(x+2)-*(xc+2)));

    *(a+0) += a4x;
    *(a+1) += a4y;
    *(a+2) += a4z;

}
*/
/* ---------- advancement ---------- */
void do_step(double dt, double *x, double *v, double *parameter) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;

	hh = dt*0.5;

	getforce(x,v,acc0, parameter);
	for (k=0;k<3;k++) {                /* first half-step */
		xt1[k] = *(x+k)+hh**(v+k);
		vt1[k] = *(v+k)+hh**(acc0+k);
	}

	getforce(&xt1[0], &vt1[0], acc1, parameter);
	for (k=0;k<3;k++) {                /* second half-step */
		xt2[k] = *(x+k)+hh*vt1[k];
		vt2[k] = *(v+k)+hh**(acc1+k);
	}

	getforce(&xt2[0], &vt2[0], acc2, parameter);
	for (k=0;k<3;k++) {                /* third half-step with results of second half-step */
		xt3[k] = *(x+k)+dt*vt2[k];
		vt3[k] = *(v+k)+dt**(acc2+k);
	}

	getforce(&xt3[0], &vt3[0], acc3, parameter);
	for (k=0;k<3;k++) {                /* Runge-Kutta formula */
		*(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
		*(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
	}

}
