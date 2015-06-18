#include "orbit.h"
#include "orbit_utils.h"
#include <Python.h>
#include <numpy/arrayobject.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int orbit(PyDictObject *input_parameters){
    //round simulation time to multiple of output time
    struct Gal gal[ngals];
    struct Params parameters;

    parameters.b1_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "b1_LMJ"));
    parameters.M1_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "M1_LMJ"));
    parameters.a2_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "a1_LMJ"));
    parameters.b2_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "b2_LMJ"));
    parameters.M2_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "M2_LMJ"));
    parameters.Mhalo = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "Mhalo"));
    parameters.q_halo = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "q_halo"));
    parameters.r_halo = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "r_halo"));

    int ratio, n;
    double tpast = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "tpast"));
    printf("%10.5f\n", tpast);
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;

    //convert coordinates into Cartesian frame
    PyArrayObject *distance_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "distance_gal");
    PyArrayObject *pm_mu_alphacosdelta = (PyArrayObject*)PyDict_GetItemString(input_parameters, "pm_mu_alphacosdelta");
    PyArrayObject *pm_mu_delta = (PyArrayObject*)PyDict_GetItemString(input_parameters, "pm_mu_delta");
    PyArrayObject *l = (PyArrayObject*)PyDict_GetItemString(input_parameters, "l");
    PyArrayObject *b = (PyArrayObject*)PyDict_GetItemString(input_parameters, "b");
    PyArrayObject *mass_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "mass_gal");
    PyArrayObject *rad_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "rad_gal");

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
    for (n=0; n<ngals; n++){
        dsuntemp = *(double*)PyArray_GETPTR1(distance_gal,n)*1000.;
        mu_alphatemp = 0.0;
        mu_alphacosdeltatemp = *(double*)PyArray_GETPTR1(pm_mu_alphacosdelta,n)*1000.0;
        mu_deltatemp = *(double*)PyArray_GETPTR1(pm_mu_delta, n)*1000.0;
        vrsuntemp = vLSR;
        ltemp = *(double*)PyArray_GETPTR1(l, n);
        btemp = *(double*)PyArray_GETPTR1(b, n); 
        lcosbtemp = 0.0;
        vLSRtemp = vLSR;
        convert(gal[n].pos,
                gal[n].vel,
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
        gal[n].mhalo = *(double*)PyArray_GETPTR1(mass_gal, n);
        gal[n].r_halo = *(double*)PyArray_GETPTR1(rad_gal, n);
    }
    //get position of cluster at t = -tpast
    double sign = -1.0;
    double tmax = tpast;
    double dtoutt = tpast;
    double t = tstart;
    int err;
    err = rk4_drv(&t, tmax, dtoutt, mdiff, gal, parameters, sign);

    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
    sign = 1.0;
    dtoutt = dtout;
    tmax = tfuture;
    err = rk4_drv(&t, tmax, dtoutt, mdiff, gal, parameters, sign);
    printf("%d", err);
    //}

	return 0;
}

/* --------------- extrapolation method --------------- */
int rk4_drv(double *t,
            double tmax,
            double dtout,
            double mdiff,
            struct Gal *gal,
            struct Params parameters,
            double sign){

	double tout, diff, dt = 0.0;
	double xe1[3], ve1[3], difftemp;
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
                        (*(gal+n)).post[k] = (*(gal+n)).pos[k];
                        (*(gal+n)).velt[k] = (*(gal+n)).vel[k];
			xe1[k]=(*(gal+n)).pos[k];
			ve1[k]=(*(gal+n)).vel[k];
		    }
		    do_step(dt, xe1, ve1, n, gal, parameters);      /* One full step */

		    do_step(0.5*dt, (*(gal+n)).post, (*(gal+n)).velt, n, gal, parameters);  /* Two half steps */
		    do_step(0.5*dt, (*(gal+n)).post, (*(gal+n)).velt, n, gal, parameters);
		    difftemp = sqrt(pow(xe1[0] - (*(gal+n)).post[0],2) +
		                    pow(xe1[1] - (*(gal+n)).post[1],2) +
		                    pow(xe1[2] - (*(gal+n)).post[2],2));
                    if (difftemp > diff) {diff = difftemp;} // hold highest value to compare with mdiff below
                }
                // end loop over each galaxy
		if (diff<=mdiff) {         /* Is difference below accuracy threshold? */
		    *t+=dt;
		    //update the test particles here
                    for (n=0; n<ngals; n++){
		        for (k=0;k<3;k++) {
			    (*(gal+n)).pos[k]=(*(gal+n)).post[k];          /* If yes -> continue and double step size */
			    (*(gal+n)).vel[k]=(*(gal+n)).velt[k];
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
	    for (n=0; n<ngals; n++){
                if ((sqrt(pow((*(gal+n)).pos[0],2)+pow((*(gal+n)).pos[1],2)+pow((*(gal+n)).pos[2],2)) < Rgalmin) ||
                    (sqrt(pow((*(gal+n)).pos[0],2)+pow((*(gal+n)).pos[1],2)+pow((*(gal+n)).pos[2],2)) > Rgalmax)  ) return 1;  /* Abort if r is too extreme */
            }
	    if (sign**t>=sign*(tout)) {
                // save info on both gals
                for (n=0; n<ngals; n++){
                    printf("Galaxy %d: x,y,z, t= %lf, %lf, %lf %lf\n", n, (*(gal+n)).pos[0], (*(gal+n)).pos[1], (*(gal+n)).pos[2], *t);
                }
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
void getforce(double *x, double *v, double *a, struct Params parameters){
	double r1, r2, r3;
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        
    // Potential for MW
	//Hernquist bulge
	r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));

	//ax += -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**x/r1;
	//ay += -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+1)/r1;
	//az += -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+2)/r1;
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
	az += -G*parameters.M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(pow(*(x+2), 2) + pow(parameters.b2_LMJ, 2)))/
	        sqrt(pow(*(x+2), 2) + pow(parameters.b2_LMJ, 2)) * *(x+2);

	//NFW Halo
	r3 = sqrt(pow(*x, 2) + pow(*(x+1), 2) + pow(*(x+2)/q_halo, 2));

        ax += -G*parameters.Mhalo/r3 * (log(1.0 + r3/parameters.r_halo)/r3 - 1.0/(parameters.r_halo+r3)) * *x/r3;
        ay += -G*parameters.Mhalo/r3 * (log(1.0 + r3/parameters.r_halo)/r3 - 1.0/(parameters.r_halo+r3)) * *(x+1)/r3;
        az += -G*parameters.Mhalo/r3 * (log(1.0 + r3/parameters.r_halo)/r3 - 1.0/(parameters.r_halo+r3)) * *(x+2)/(parameters.q_halo*parameters.q_halo*r3);

	*(a+0) = ax;
	*(a+1) = ay;
	*(a+2) = az;
}


void getforce_gals(double *x, double *v, double *a, int gal_num, struct Gal *gal, struct Params parameters){
    getforce(x, v, a, parameters);
    int i;
    double r;
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
//    double r200, c, deltachar;
//    double k = 1.3e-7;
    //NFW Halo
    for (i=0; i<ngals; i++){
        if (i != gal_num){
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
                     pow(*(x+1) - gal[i].pos[1], 2) +
                     pow(*(x+2) - gal[i].pos[2], 2));
        
            ax += -G* gal[i].mhalo/r * (log(1.0 + r/ gal[i].r_halo)/r - 1.0/( gal[i].r_halo+r)) * (*x-gal[i].pos[0])/r;
            ay += -G* gal[i].mhalo/r * (log(1.0 + r/ gal[i].r_halo)/r - 1.0/( gal[i].r_halo+r)) * (*(x+1)-gal[i].pos[1])/r;
            az += -G* gal[i].mhalo/r * (log(1.0 + r/ gal[i].r_halo)/r - 1.0/( gal[i].r_halo+r)) * (*(x+2)-gal[i].pos[2])/r;
        }
    }
    *(a+0) += ax;
    *(a+1) += ay;
    *(a+2) += az;

}


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
void do_step(double dt, double *x, double *v, int gal_num, struct Gal *gal, struct Params parameters) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;

	hh = dt*0.5;
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
