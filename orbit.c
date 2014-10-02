#include "orbit.h"
#include "orbit_utils.h"
#include "orbit_data.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int orbit(double *fit,
          double mass_cluster,
          double pm_mu_delta,
          double pm_mu_alphacosdelta,
          double mass_halo,
          double distance_cluster,
          double mass_loss_rate,
          double q_halo,
          double r_halo,
          double tpast,
          double rgalsun,
          double vLSR,
          double sigma_x,
          double sigma_v,
          double sigma_vx,
          double sigma_mu,
          int newspaper) {


	*fit = SUPERSMALL;
	int err;

    //round simulation time to multiple of output time
    int ratio;
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;

    //put the MCMC parameters in the appropriate array for passing it to the integrator
    double parameter[17];
    parameter[0] = mass_cluster;            //mass of cluster
    parameter[1] = pm_mu_delta;             //delta component of proper motion
    parameter[2] = pm_mu_alphacosdelta;     //alpha component of proper motion
    parameter[3] = mass_halo;               //mass of DM halo
    parameter[4] = distance_cluster;        //distance cluster-sun
    parameter[6] = rgalsun;                 //distance sun-galactic center
    parameter[7] = vLSR;                    //y-velocity of LSR
    parameter[8] = mass_loss_rate;               //cluster mass loss rate
    parameter[10] = q_halo;                 //flattening of dark halo
    parameter[11] = r_halo;                 //concentration of dark halo
    parameter[12] = tpast;                //integration time
    parameter[13] = sigma_x;                //smoothing position of overdensities
    parameter[14] = sigma_v;                //smoothing v_r
    parameter[15] = sigma_vx;               //smoothing position of v_r stars
    parameter[16] = sigma_mu;              //smoothing proper motions of v_r stars


    //convert coordinates into Cartesian frame
	double x[3];
    double v[3];
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
	dsuntemp = parameter[4]*1000.0;
	mu_alphatemp = 0.0;
	mu_alphacosdeltatemp = parameter[2]*1000.0;
	mu_deltatemp = parameter[1]*1000.0;
	vrsuntemp = vr;
	ltemp = l;
	btemp = b;
	lcosbtemp = lcosb;
	vLSRtemp = vLSR;
	convert(x,
            v,
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
	getforce(x,v,atemp,parameter);


	//get the acceleration at solar circle
	double xsuntemp[3], vsuntemp[3], asuntemp[3], Vcirc;
	xsuntemp[0]=rgalsun;
	xsuntemp[1]=0.0;
	xsuntemp[2]=0.0;
	vsuntemp[0]=0.0;
	vsuntemp[1]=0.0;
	vsuntemp[2]=0.0;
	getforce(xsuntemp,vsuntemp,asuntemp,parameter);
	Vcirc = sqrt(rgalsun*sqrt(asuntemp[0]*asuntemp[0]+asuntemp[1]*asuntemp[1]+asuntemp[2]*asuntemp[2]));

    //set unrealistic parameters to artificially low likelihood value
	//if (Vcirc > 300.0 || Vcirc < 160.0) {
    //    *fit = SUPERSMALL;
    //} else {

    //get position of cluster at t = -tpast
    double sign = -1.0;
    double tmax = tpast;
    double dtoutt = tpast;
    double t = tstart;
    rk4_drv(&t,tmax,dtoutt,mdiff,&x[0],&v[0],sign,parameter,fit,newspaper);


    //save whatever initial conditions for possible restarts
    double xrestart[3], vrestart[3], trestart;
    xrestart[0] = x[0];
    xrestart[1] = x[1];
    xrestart[2] = x[2];
    vrestart[0] = v[0];
    vrestart[1] = v[1];
    vrestart[2] = v[2];
    trestart = t;

    double tspan = sqrt(pow(tpast-tstart,2));
    *fit = SUPERSMALL;

    //integrate cluster orbit backwards to t = -tint
    x[0] = xrestart[0];
    x[1] = xrestart[1];
    x[2] = xrestart[2];
    v[0] = vrestart[0];
    v[1] = vrestart[1];
    v[2] = vrestart[2];
    t = trestart;

    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
    sign = 1.0;
    dtoutt = dtout;
    tmax = tfuture;
    if (radio) printf("\nt = %f\tdtout = %f\n",t,dtout);
    err = rk4_drv(&t,tmax,dtoutt,mdiff,&x[0],&v[0],sign,parameter,fit,newspaper);

    //}

    //set lower bound on likelihood value
	if (err || *fit < SUPERSMALL) *fit = SUPERSMALL;

	printf("%9.4f\t%7.1f\t%9.7f\t%9.7f\t%10.6e\t%5.3f\t%8.1f\t%6.4f\t%6.1f\t%10.4f\t%6.1f\t%6.1f\t%6.2f\t%.0f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n",
            *fit, parameter[0], parameter[1], parameter[2], parameter[3],
            parameter[10], parameter[11], parameter[4], parameter[8], parameter[6],
            parameter[7], sqrt(atemp[0]*atemp[0]+atemp[1]*atemp[1]+atemp[2]*atemp[2]),
            Vcirc, tpast, sigma_x, sigma_v, sigma_vx, sigma_mu);

	return 0;
}


/* --------------- extrapolation method --------------- */
int rk4_drv(double *t,
            double tmax,
            double dtout,
            double mdiff,
            double *x,
            double *v,
            double sign,
            double *parameter,
            double *fit,
            int newspaper){

    void write_snapshot_tail(double next_snap, double *xtemp, double *vtemp);
    void write_snapshot_cluster(double next_snap, double *xtemp, double *vtemp);
	double tout,diff,dt, mdifft,dttemp = 0.0;
    double next_snap_cluster, next_snap_tail = 0.0;
	double atemp[3], xe1[3], xe2[3], ve1[3], ve2[3], xt[3], vt[3], vmaxwell,vtemp, omega[3], omegat, xc[3], vc[3], ac[3];
	double actualclustermass;
    double rgalsun = parameter[6];
	int k,i;
	int err;
	double tt, r;
	mdifft = mdiff;
    char name[50];
    FILE *fz;
    sprintf(name, "temp.dat");
    double tpast;
    tpast = parameter[12];
    if (newspaper) fz = fopen(name,"w");
    //int snapnumber = 0;

	//initialize output/insertion of tail particles
	tout = *t;							/* time of next output/insertion */
    next_snap_cluster = *t;
    next_snap_tail = *t;
	dt = sign*dt0;                /* initial time step */
    int NMAX = 2.0*(tstart-tpast)/(1.0*dtout);

	int columns = 5;
	double **star;
	star = (double **)calloc(NMAX,sizeof(double *));
	for (i=0;i<NMAX;i++){
		star[i] = (double *)calloc(columns,sizeof(double));
		if (star[i] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}

	int starcount = 0;
	*fit = SUPERSMALL;


	//integrate cluster
	do {

		/***********
		 * CLUSTER *
		 ***********/

		if (sign**t>=sign*(tout)) {              /* output */

			//coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
			//vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities
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
			vLSRtemp = parameter[7];


			//integrate test tail particle for leading and trailing tail from current cluster position to time of observation, i.e. today
			if ((tails) && (sign>0.0) && (sign**t<=sign*tmax) && (sign**t<tstart))  {

				actualclustermass = parameter[0]+sqrt(*t**t)*parameter[8];
				parameter[9] = sqrt(*t**t);

				if (radio) printf("M(t) = %f\n", actualclustermass);


				/****************
				 * LEADING TAIL *
				 ****************/

				r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
				vtemp = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
				omega[0] = (x[1]*v[2]-x[2]*v[1])/(r*r);
				omega[1] = (x[2]*v[0]-x[0]*v[2])/(r*r);
				omega[2] = (x[0]*v[1]-x[1]*v[0])/(r*r);
				omegat = sqrt(omega[0]*omega[0]+omega[1]*omega[1]+omega[2]*omega[2]);

				double rtide, rtidef, dphi;

				double at[3], at2[3]; //force evaluation for dphi
				for (i=0;i<3;i++) {
					xt[i] = x[i]/(1.0*r)*(r-20.0);
					vt[i] = 0.0;
				}
				getforce(xt,vt,at,parameter);
				for (i=0;i<3;i++) xt[i] = x[i]/(1.0*r)*(r+20.0);
				getforce(xt,vt,at2,parameter);

				dphi = (sqrt(at[0]*at[0]+at[1]*at[1]+at[2]*at[2])-sqrt(at2[0]*at2[0]+at2[1]*at2[1]+at2[2]*at2[2]))/40.0;
				rtide = pow(G*parameter[0]/sqrt(pow((dphi+omegat*omegat),2)),1.0/3.0);


				if (rtide>parameter[5]) rtidef = rtide; //in case edge radius is wanted
				else rtidef = parameter[5];

				do {
					tt = *t;
					for (i=0;i<3;i++) {
						xc[i] = x[i];
						vc[i] = v[i];
					}

					if (radio) printf("tt = %.0f\trtide = %.1f\tomega = %.5f\tr = %.1f\t\t", tt, rtidef,omegat,r);

					vmaxwell = 1.0/3.0*sqrt(pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2));
					//printf("\tv1 = %g", vmaxwell);

					for (i=0;i<3;i++) {
                        xt[i] = x[i]/(1.0*r)*(r-rtidef)+rexcess*rtidef*get_gauss();
                        vt[i] = v[i]/vtemp*(vtemp-omegat*rtidef)-(vmaxwell*x[i]/(1.0*r));
					}

					err = rk4_drvtail(&tt,tstart,mdifft,&xt[0],&vt[0],sign,xc,vc,parameter, next_snap_tail);

					if (err==1) {
						for (i=0;i<NMAX;i++) free (star[i]);
						free(star);
                        if (newspaper) fclose(fz);
						return 1;
					} else if (err ==2) {
						rtidef += 10.0;
						if (rtidef>rtidemax) {
							for (i=0;i<NMAX;i++) free (star[i]);
							free(star);
							return 1;
						}
					}
				} while (err);
                //if (snapshot) write_snapshot(snapfile, xt, vt);

                //store galactic coordinates and radial velocity of tail particle
				//coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
				//vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities
				convert(xt, vt, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp,
                        &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp,rgalsun);

                if (newspaper) fprintf(fz,"%lf %lf %lf %lf %lf %lf\n", ltemp, btemp, lcosbtemp, vrsuntemp, RAtemp, DECtemp);

				star[starcount][0] = RAtemp;
				star[starcount][1] = DECtemp;
				star[starcount][2] = vrsuntemp;
                star[starcount][3] = mu_alphacosdeltatemp;
                star[starcount][4] = mu_deltatemp;
				starcount++;
                next_snap_tail = *t;




				/*****************
				 * TRAILING TAIL *
				 *****************/

				do {
					tt = *t;
					for (i=0;i<3;i++) {
						xc[i] = x[i];
						vc[i] = v[i];
					}

					vmaxwell = 1.0/3.0*sqrt(pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2));
					//printf("\tv2 = %g\n", vmaxwell);

					for (i=0;i<3;i++) {
                        xt[i] = x[i]/r*(r+rtidef)+rexcess*rtidef*get_gauss();
                        vt[i] = v[i]/vtemp*(vtemp+omegat*rtidef)+(vmaxwell*x[i]/(1.0*r));
					}

					err = rk4_drvtail(&tt,tstart,mdifft,&xt[0],&vt[0],sign,xc,vc,parameter, next_snap_tail);

					if (err==1) {
						for (i=0;i<NMAX;i++) free (star[i]);
						free(star);
                        if (newspaper) fclose(fz);
						return 1;
					} else if (err ==2) {
						rtidef += 10.0;
						if (rtidef>rtidemax) {
							for (i=0;i<NMAX;i++) free (star[i]);
							free(star);
							return 1;
						}
					}
				} while (err);

                //store galactic coordinates and radial velocity of tail particle
				//coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
				//vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities
				//if (snapshot) write_snapshot(snapfile, xt, vt);
                convert(xt, vt, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp,
                        &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp, rgalsun);

                if (newspaper) fprintf(fz,"%lf %lf %lf %lf %lf %lf\n", ltemp, btemp, lcosbtemp, vrsuntemp, RAtemp, DECtemp);

				star[starcount][0] = RAtemp;
				star[starcount][1] = DECtemp;
				star[starcount][2] = vrsuntemp;
                star[starcount][3] = mu_alphacosdeltatemp;
                star[starcount][4] = mu_deltatemp;
				starcount++;
                next_snap_tail = *t;
			}


			tout+=dtout;              /* increase time of output/next insertion */
		}



        //advance cluster particle
		int count = 0;
        int laststep = 0;
		do {
			if (sign*(*t+dt) > sign*tout) {
				dt = tout-*t;
                laststep = 1;
			}
			for (k=0;k<3;k++) {
				xe1[k]=x[k];
				xe2[k]=x[k];
				ve1[k]=v[k];
				ve2[k]=v[k];
			}
			do_step(dt,xe1,ve1,parameter);      /* One full step */

			do_step(0.5*dt,xe2,ve2,parameter);  /* Two half steps */
			do_step(0.5*dt,xe2,ve2,parameter);

			diff = sqrt(pow(*xe1 - *xe2,2) + pow(*(xe1+1) - *(xe2+1),2) + pow(*(xe1+2) - *(xe2+2),2));

			if (diff<mdiff) {         /* Is difference below accuracy threshold? */
				*t+=dt;
				dttemp = dt;

				for (k=0;k<3;k++) {
					x[k]=xe2[k];          /* If yes -> continue and double step size */
					v[k]=ve2[k];
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
        if (snapshot)
        {
            if (sign**t>=sign*(next_snap_cluster))
                {
                    write_snapshot_cluster(next_snap_cluster, xe2, ve2);
                    next_snap_cluster += dtsnap;
                }
		}
        if ((sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) < Rgalmin) || (sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) > Rgalmax)  ) return 1;  /* Abort if r is too extreme */
	} while (sign**t<sign*(tmax));


	//last printout
	double dsuntemp, vrsuntemp, vrtemp, ltemp, btemp, lcosbtemp, RAtemp, DECtemp, mu_alphatemp, mu_alphacosdeltatemp, mu_deltatemp, mutemp, PAtemp, vLSRtemp;
	vLSRtemp = parameter[7];
	getforce(x,v,atemp,parameter);
	convert(x, v, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp,rgalsun);
	if (radio) printf("\n%.2f\n\n", sqrt(atemp[0]*atemp[0]+atemp[1]*atemp[1]+atemp[2]*atemp[2]));

    if (newspaper) fclose(fz);
    //if (snapshot) fclose(snapfile);
    if (starcount) *fit = fitness(starcount, star, parameter[13], parameter[14], parameter[15], parameter[16]);

	for (i=0;i<NMAX;i++) free (star[i]);
	free(star);

	return 0;

}

/* --------------- Write Snapshot (currently in csv format) --------------- */
void write_snapshot_tail(double next_snap, double *xtemp, double *vtemp){
    double x = xtemp[0] / 1000.0;
    double y = xtemp[1] / 1000.0; // keep things in kpc
    double z = xtemp[2] / 1000.0;
    double v = sqrt(vtemp[0] * vtemp[0] + vtemp[1] * vtemp[1] + vtemp[2] * vtemp[2]);
    FILE *snapfile;
    char snapname[50];
    int snapnumber = - next_snap / dtsnap; // integer number of snapshots
    struct stat st;
    sprintf(snapname, "snapshot_tail.csv.%i", snapnumber);
    int result = stat(snapname, &st);
    if (!result)
    {
       snapfile = fopen(snapname,"a");
    } else {
        snapfile = fopen(snapname,"w");
        fprintf(snapfile,"X,Y,Z,V\n");
    }
    fprintf(snapfile,"%10.5f,%10.5f,%10.5f,%10.5f\n", x, y, z, v);
    fclose(snapfile);
}

/* --------------- Write Snapshot (currently in csv format) --------------- */
void write_snapshot_cluster(double next_snap, double *xtemp, double *vtemp){
    double x = xtemp[0] / 1000.0;
    double y = xtemp[1] / 1000.0;
    double z = xtemp[2] / 1000.0;
    double v = sqrt(vtemp[0] * vtemp[0] + vtemp[1] * vtemp[1] + vtemp[2] * vtemp[2]);
    FILE *snapfile;
    char snapname[50];
    int snapnumber = next_snap / dtsnap; // integer number of snapshots
    sprintf(snapname, "snapshot_cluster.csv.%i", snapnumber);
    snapfile = fopen(snapname,"w");
    fprintf(snapfile,"X,Y,Z,V\n");
    fprintf(snapfile,"%10.5f,%10.5f,%10.5f,%10.5f\n", x, y, z, v);
    fclose(snapfile);

}


/* --------------- extrapolation method tail --------------- */
int rk4_drvtail(double *t, double tmax, double mdiff, double *x, double *v, double sign, double *xc, double *vc, double *parameter, double next_snap_tail){
	double diff,dt;
	double xe1[3], xe2[3], ve1[3], ve2[3];
	double xce1[3], xce2[3], vce1[3],vce2[3];
	int k;
    int laststep = 0;

	dt = sign*dt0;		/* initial time step */

	while (sign**t<sign*tmax) {
		if (sign*(*t+dt) > sign*tmax) {
			dt = tmax-*t;
            laststep = 1;
		}

		do {
			for (k=0;k<3;k++) {
				xe1[k]=x[k];
				xe2[k]=x[k];
				ve1[k]=v[k];
				ve2[k]=v[k];
				xce1[k]=xc[k];
				xce2[k]=xc[k];
				vce1[k]=vc[k];
				vce2[k]=vc[k];
			}

			parameter[9] = sqrt(*t * *t);

			do_steptail(dt,xe1,ve1,xce1,vce1, parameter);

			do_steptail(0.5*dt,xe2,ve2,xce2,vce2, parameter);
			do_steptail(0.5*dt,xe2,ve2,xce2,vce2, parameter);

			diff = sqrt(pow(*xe1 - *xe2,2) + pow(*(xe1+1) - *(xe2+1),2) + pow(*(xe1+2) - *(xe2+2),2));

			if (diff<mdiff) {
				*t+=dt;

				for (k=0;k<3;k++) {
					x[k]=xe2[k];
					v[k]=ve2[k];
					xc[k]=xce2[k];
					vc[k]=vce2[k];
				}
				dt = dt*2.0;

			} else {
				dt = dt/2.0;
			}

			if ((sign*dt < 0.01*dt0) && !laststep) {
				printf("Aborted... dt = %lf\n", dt);
				*t = tmax*2.0;
				diff = mdiff/2;
				return 1;
			}

		} while (diff>mdiff);
        if (snapshot)
        {
            if (sign**t>=sign*(next_snap_tail))
                {
                    write_snapshot_tail(next_snap_tail, xe2, ve2);
                    next_snap_tail += dtsnap;
                }
        }

		if (sqrt(pow(x[0]-xc[0],2)+pow(x[1]-xc[1],2)+pow(x[2]-xc[2],2)) < Rstop) return 2; //increase Rtide
		else if (sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) < Rgalmin) return 1;
        //Add snapshot here??
	}

	return 0;

}

/* ----------- force ----------- */
void getforce(double *x, double *v, double *a, double *parameter){
	double r1, r2, r3;
	double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
    double q1, q2, qz, phi, C1, C2, C3, rhalo, vhalo;
//    double r200, c, deltachar;
//    double k = 1.3e-7;

	if (gpot == 1) {
		//Allen & Santillan (1991) potential w updated values from Irrgang et al. (2013)

		//Point mass
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2) + b1*b1);

		a1x = -G*M1/(r1*r1*r1)**x;
		a1y = -G*M1/(r1*r1*r1)**(x+1);
		a1z = -G*M1/(r1*r1*r1)**(x+2);

		//Miyamato disk
		r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2 + sqrt(*(x+2) * *(x+2) + b2*b2),2));

		a2x = -G*M2/(r2*r2*r2) * *x;
		a2y = -G*M2/(r2*r2*r2) * *(x+1);
		a2z = -G*M2/(r2*r2*r2) * (a2 + sqrt(*(x+2) * *(x+2) + b2*b2))/sqrt(*(x+2) * *(x+2) + b2*b2) * *(x+2);

		//Log Halo
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);

		a3x = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *x/r3;
		a3y = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *(x+1)/r3;
		a3z = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *(x+2)/(parameter[10]*parameter[10]*r3);

	} else if (gpot == 2) {
		//Log Halo from Koposov et al. (2010)
		r3 = *x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]; //R^2!

		a3x = -G*parameter[3]/parameter[11] *  *x/r3;
		a3y = -G*parameter[3]/parameter[11] *  *(x+1)/r3;
		a3z = -G*parameter[3]/parameter[11] *  *(x+2)/(parameter[10]*parameter[10]*r3);

		a1x = a1y = a1z = a2x = a2y = a2z = 0.0;

	} else if (gpot == 3) {
		//potential from Johnston/Law/Majewski/Helmi

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
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);

        a3x = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *x/r3;
        a3y = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+1)/r3;
        a3z = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+2)/(parameter[10]*parameter[10]*r3);

//        a3x = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *x/r3;
//		a3y = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+1)/r3;
//		a3z = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+2)/(parameter[10]*parameter[10]*r3);

	} else if (gpot == 4) {
        //potential from Law+Majewski 2010

        //Hernquist bulge
        r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));

        a1x = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**x/r1;
        a1y = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+1)/r1;
        a1z = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+2)/r1;

        //Miyamato disk
        r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ),2));

        a2x = -G*M2_LMJ/(r2*r2*r2) * *x;
        a2y = -G*M2_LMJ/(r2*r2*r2) * *(x+1);
        a2z = -G*M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))/sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ) * *(x+2);

        //Logarithmic Halo
        //USE FIXED PARAMETERS FROM PEARSON+2014, LM+2010
        q2 = 1.0;
        // Pay attention!!
        vhalo = parameter[3];  // this is the "MHALO" Parameter for all other halos... for this halo it is the Vhalo!!
        rhalo = parameter[11];
        phi = 97*PI/180.;
        q1 = 1.38;
        qz = parameter[10];  // "q_halo"
        C1 = pow(cos(phi), 2)/pow(q1, 2) + pow(sin(phi), 2)/pow(q2, 2);
        C2 = pow(cos(phi), 2)/pow(q2, 2) + pow(sin(phi), 2)/pow(q1, 2);
        C3 = 2 * sin(phi) * cos(phi) * (1.0/pow(q1, 2) - 1.0/pow(q2, 2));
        double vh2 = vhalo*vhalo;

        r3 = C1* *x * *x + C2* *(x+1) * *(x+1) + C3* *x * *(x+1) + *(x+2)/qz * *(x+2)/qz;

        double fac = 0.5*vh2/(r3+rhalo*rhalo);
        a3x = -fac*(2.*C1* *x +C3 * *(x+1));
        a3y = -fac*(2.*C2* *(x+1)+C3* *x);
        a3z = -fac*(2.* *(x+2)/qz/qz);

    }

	*(a+0) = a1x + a2x + a3x;
	*(a+1) = a1y + a2y + a3y;
	*(a+2) = a1z + a2z + a3z;

}

/* ----------- force tail ----------- */
void getforcetail(double *x, double *v, double *a, double *xc, double *parameter){
	//double r1, r2, r3, 
    double r4;
	//double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, 
    double a4x, a4y, a4z;
	double actualmass;
    //double q1, q2, qz, phi, C1, C2, C3, rhalo, vhalo;
	actualmass = parameter[0]+parameter[9]*parameter[8];

    getforce(x, v, a, parameter);

	//Cluster
	r4 = sqrt((*x-*xc) * (*x-*xc) + (*(x+1)-*(xc+1)) * (*(x+1)-*(xc+1)) + (*(x+2)-*(xc+2)) * (*(x+2)-*(xc+2)) + R4*R4);

	a4x = -G*actualmass/(r4*r4*r4)*(*x-*xc);
	a4y = -G*actualmass/(r4*r4*r4)*(*(x+1)-*(xc+1));
	a4z = -G*actualmass/(r4*r4*r4)*(*(x+2)-*(xc+2));
/*
	*(a+0) = a1x + a2x + a3x + a4x;
	*(a+1) = a1y + a2y + a3y + a4y;
	*(a+2) = a1z + a2z + a3z + a4z;
*/

    *(a+0) += a4x;
    *(a+1) += a4y;
    *(a+2) += a4z;

}

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

/* ---------- advancement tail ---------- */
void do_steptail(double dt, double *x, double *v, double *xc, double *vc, double *parameter) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3] ,xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	double ac0[3],ac1[3],ac2[3],ac3[3],xct1[3],xct2[3],xct3[3],vct1[3],vct2[3],vct3[3];
	int k;

	hh = dt*0.5;

	getforcetail(x,v,acc0,xc,parameter);
	getforce(xc,vc,ac0,parameter);
	for (k=0;k<3;k++) {
		xt1[k] = *(x+k)+hh**(v+k);
		vt1[k] = *(v+k)+hh**(acc0+k);
		xct1[k] = *(xc+k)+hh**(vc+k);
		vct1[k] = *(vc+k)+hh**(ac0+k);
	}

	getforcetail(&xt1[0], &vt1[0], acc1, xct1, parameter);
	getforce(&xct1[0], &vct1[0], ac1, parameter);
	for (k=0;k<3;k++) {
		xt2[k] = *(x+k)+hh*vt1[k];
		vt2[k] = *(v+k)+hh**(acc1+k);
		xct2[k] = *(xc+k)+hh*vct1[k];
		vct2[k] = *(vc+k)+hh**(ac1+k);
	}

	getforcetail(&xt2[0], &vt2[0], acc2, xct2, parameter);
	getforce(&xct2[0], &vct2[0], ac2, parameter);
	for (k=0;k<3;k++) {
		xt3[k] = *(x+k)+dt*vt2[k];
		vt3[k] = *(v+k)+dt**(acc2+k);
		xct3[k] = *(xc+k)+dt*vct2[k];
		vct3[k] = *(vc+k)+dt**(ac2+k);
	}

	getforcetail(&xt3[0], &vt3[0], acc3, xct3, parameter);
	getforce(&xct3[0], &vct3[0], ac3, parameter);
	for (k=0;k<3;k++) {
		*(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
		*(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
		*(xc+k) += dt/6.0*(*(vc+k)+2.0*(vct1[k]+vct2[k])+vct3[k]);
		*(vc+k) += dt/6.0*(*(ac0+k)+2.0*(*(ac1+k)+*(ac2+k))+*(ac3+k));
	}

}

/* ---------- likelihood evaluation ----------- */
double fitness(int N, double **star, double sigma_x, double sigma_v, double sigma_vx, double sigma_mu){
	if(radio) printf("\nGetting likelihood value...\n");
	int i, j;
    double chi2;
    double sigma_x2, sigma_v2, sigma_vx2, sigma_mu2;
    double n_OD_NGC5466(int i, int j);
    double n_vr_NGC5466(int i, int j);
    double n_Belo_NGC5466(int i, int j);
    double n_Lux_NGC5466(int i, int j);
    sigma_x2 = sigma_x*sigma_x;
    sigma_v2 = sigma_v*sigma_v;
    sigma_vx2 = sigma_vx*sigma_vx;
    sigma_mu2 = sigma_mu*sigma_mu;


    int nr_OD;
    if (!data_type) {
        nr_OD = 20; //   3 by hand from tails of Belokurov + 17 from Lux
    } else {
        nr_OD = 3;    // 3 from Belokurov
    }

    double n_OD[nr_OD][6];

    if (!data_type) {
        for (i=0;i<nr_OD;i++){
            for (j=0;j<6;j++){
                n_OD[i][j] = n_OD_NGC5466(i, j);
            }
        }
    } else {
        // Lux points
        for (i=0;i<nr_OD;i++){
            for (j=0;j<6;j++){
                n_OD[i][j] = n_Belo_NGC5466(i, j);
            }
        }
    }

    int nr_vr = 6;
    double n_vr[nr_vr][10];

    for (i=0;i<nr_vr;i++){
        for (j=0;j<10;j++){
            n_vr[i][j] = n_vr_NGC5466(i, j);
        }
    }

	double dl, db, dl2, db2, dvx, dvx2;
	dl = sqrt(0.125*0.125+sigma_x2);
	db = sqrt(0.125*0.125+sigma_x2);
	dvx = sqrt(0.125*0.125+sigma_vx2);
	dl2 = dl*dl; //squared values
	db2 = db*db;
    dvx2 = dvx*dvx;

    for (j=0;j<nr_OD;j++) {
        n_OD[j][5] =  sigma_x2+n_OD[j][4]*n_OD[j][4]/(8.0*log(2.0));//squared values of sigma_x plus sigma_obs
    }

	for(i=0;i<N;i++) {

		//if (star[i][0]>180.0) star[i][0] = star[i][0]-360.0; for use with galactic coordinates

        //velocities
        for (j=0;j<nr_vr;j++) {
            n_vr[j][5] += exp(-0.5*((star[i][0]-n_vr[j][4])*(star[i][0]-n_vr[j][4])/dvx2 +
                                    (star[i][1]-n_vr[j][3])*(star[i][1]-n_vr[j][3])/dvx2 +
                                    (star[i][2]-n_vr[j][0])*(star[i][2]-n_vr[j][0])/(n_vr[j][1]*n_vr[j][1]+sigma_v2) +
                                    (star[i][3]-n_vr[j][6])*(star[i][3]-n_vr[j][6])/(n_vr[j][8]*n_vr[j][8]+sigma_mu2) +
                                    (star[i][4]-n_vr[j][7])*(star[i][4]-n_vr[j][7])/(n_vr[j][9]*n_vr[j][9]+sigma_mu2) ));
        }

        //overdensities
        for (j=0;j<nr_OD;j++) {
            n_OD[j][2] +=  exp(-0.5*((star[i][0]-n_OD[j][0])*(star[i][0]-n_OD[j][0])/n_OD[j][5] +
                                (star[i][1]-n_OD[j][1])*(star[i][1]-n_OD[j][1])/n_OD[j][5]));
        }

	}

    if (N) {
        //normalization velocities
        for (j=0;j<nr_vr;j++) {
            n_vr[j][5] *= (1.0/(1.0*N*dvx2*sqrt(n_vr[j][1]*n_vr[j][1]+sigma_v2)
                                *sqrt(n_vr[j][8]*n_vr[j][8]+sigma_mu2)
                                *sqrt(n_vr[j][9]*n_vr[j][9]+sigma_mu2)));
        }

        //normalization overdensities
        for (j=0;j<nr_OD;j++) {
            n_OD[j][2] *=  (1.0/(1.0*N*n_OD[j][5]));
        }

        //construction of final likelihood value
        chi2 = 0.0;

		for (i=0;i<nr_OD;i++) chi2 += n_OD[i][3]*log(n_OD[i][2]+SMALL);
        for (i=0;i<nr_vr;i++) chi2 += log(n_vr[i][5]+SMALL);


	} else {
        chi2 = SUPERSMALL;
	}

	return chi2;

}
