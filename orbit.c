#include "orbit.h"
#include "orbit_utils.h"
#include <Python.h>
#include <numpy/arrayobject.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int orbit(int int_mode,
          int ngals,
          PyDictObject *input_parameters,
          double* output_pos,
          double* output_vel){
    // Note: no input error checking done here. Do that in python function above.    

    //Construct the parameters from the python dictionary
    struct Params parameters;
    parameters.b1_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "b1_LMJ"));
    parameters.M1_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "M1_LMJ"));
    parameters.a2_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "a1_LMJ"));
    parameters.b2_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "b2_LMJ"));
    parameters.M2_LMJ = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "M2_LMJ"));
    parameters.Mhalo = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "Mhalo"));
    parameters.q_halo = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "q_halo"));
    parameters.r_halo = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "r_halo"));
    parameters.tpast = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "tpast"));
    parameters.gamma = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "gamma_halo"));
    double dtout = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "dtout"));
    parameters.ngals = ngals;
    PyObject* folder = PyDict_GetItemString(input_parameters, "outputdir");
    char* folderstr = PyString_AsString(folder);
    parameters.outputdir = folderstr;
    // get values for our companion galaxies from the python dictionary
    PyArrayObject *mass_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "mass_gal");
    PyArrayObject *rad_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "rad_gal");
    PyArrayObject *gamma_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "gamma_gal");
    PyArrayObject *a2_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "a2_gal");
    PyArrayObject *b2_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "b2_gal");
    PyArrayObject *m2_gal = (PyArrayObject*)PyDict_GetItemString(input_parameters, "m2_gal");

    struct Gal *gal = (struct Gal *) malloc(ngals*sizeof(struct Gal));
    int ratio, n, i;
    double tpast = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "tpast"));
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;
    for (n=0; n<ngals; n++){
            gal[n].mhalo = *((double *) PyArray_GETPTR1(mass_gal, n));
            gal[n].r_halo = *((double *) PyArray_GETPTR1(rad_gal, n));
            gal[n].gamma = *((double *) PyArray_GETPTR1(gamma_gal, n));
            gal[n].a2_LMJ = *((double *) PyArray_GETPTR1(a2_gal, n));
            gal[n].b2_LMJ = *((double *) PyArray_GETPTR1(b2_gal, n));
            gal[n].M2_LMJ = *((double *) PyArray_GETPTR1(m2_gal, n));
    }
    //Set galaxies in proper place.
    PyArrayObject *pos = (PyArrayObject*)PyDict_GetItemString(input_parameters, "pos");
    PyArrayObject *vel = (PyArrayObject*)PyDict_GetItemString(input_parameters, "vel");
    for (n=0; n<ngals; n++){
        for (i=0; i<3; i++){
            gal[n].pos[i] = *(double*)PyArray_GETPTR2(pos, n, i);
            gal[n].vel[i] = *(double*)PyArray_GETPTR2(vel, n, i);
        }
    }

    //get position of cluster at t = -tpast
    if (int_mode == 1)
    {
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
    }
    else 
    {
        printf("move forward\n");
        double sign;
        if (tpast < 0.0)
            sign = -1.0;
        else 
            sign = 1.0;
        double tmax = tpast;
        double dtoutt = dtout;
        double t = tstart;
        int err;
        err = rk4_drv(&t, tmax, dtoutt, mdiff, gal, parameters, sign);
        printf("%d", err);
    }

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
            double mdiff,
            struct Gal *gal,
            struct Params parameters,
            double sign){
        int snapnum = 0;
	double tout, diff, dt = 0.0;
	double xe1[3], ve1[3], difftemp;
	int k, n;
	int ngals = parameters.ngals;
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
	
	    if (sign**t>=sign*(tout)) {
                // save info on both gals
                if (sign > 0){
                    write_snapshot(parameters, gal, *t, snapnum);
                    snapnum += 1;
                }
                tout+=dtout;              /* increase time of output/next insertion */
            }

        } while (sign**t<sign*(tmax));

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

	//Dehnen Halo
	r3 = sqrt(pow(*x, 2) + pow(*(x+1), 2) + pow(*(x+2)/q_halo, 2));
        ax += -G*parameters.Mhalo * pow(r3/(parameters.r_halo + r3), -parameters.gamma)/
                                    pow(parameters.r_halo + r3, 3) * *x;
        ay += -G*parameters.Mhalo * pow(r3/(parameters.r_halo + r3), -parameters.gamma)/
                                    pow(parameters.r_halo + r3, 3) * *(x+1);
        az += -G*parameters.Mhalo * pow(r3/(parameters.r_halo + r3), -parameters.gamma)/
                                    pow(parameters.r_halo + r3, 3) * *(x+2)/pow(parameters.q_halo, 2);
        /*
        ax += -G*parameters.Mhalo/r3 * (log(1.0 + r3/parameters.r_halo)/r3 - 1.0/(parameters.r_halo+r3)) * *x/r3;
        ay += -G*parameters.Mhalo/r3 * (log(1.0 + r3/parameters.r_halo)/r3 - 1.0/(parameters.r_halo+r3)) * *(x+1)/r3;
        az += -G*parameters.Mhalo/r3 * (log(1.0 + r3/parameters.r_halo)/r3 - 1.0/(parameters.r_halo+r3)) * *(x+2)/(parameters.q_halo*parameters.q_halo*r3);
        */
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
    int ngals = parameters.ngals;
//    double r200, c, deltachar;
//    double k = 1.3e-7;
    for (i=0; i<ngals; i++){
        if (i != gal_num){
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


            //NFW Halo
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
                     pow(*(x+1) - gal[i].pos[1], 2) +
                     pow(*(x+2) - gal[i].pos[2], 2));
            ax += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*x - gal[i].pos[0]);
            ay += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*(x+1) - gal[i].pos[1]);
            az += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*(x+2) - gal[i].pos[2]);

            
            /*
            ax += -G* gal[i].mhalo/r * (log(1.0 + r/ gal[i].r_halo)/r - 1.0/( gal[i].r_halo+r)) * (*x-gal[i].pos[0])/r;
            ay += -G* gal[i].mhalo/r * (log(1.0 + r/ gal[i].r_halo)/r - 1.0/( gal[i].r_halo+r)) * (*(x+1)-gal[i].pos[1])/r;
            az += -G* gal[i].mhalo/r * (log(1.0 + r/ gal[i].r_halo)/r - 1.0/( gal[i].r_halo+r)) * (*(x+2)-gal[i].pos[2])/r;*/

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


void write_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber){
    int n;
    int ngals = parameters.ngals;
    char *folder = parameters.outputdir;
    FILE *snapfile;
    char snapname[50];
    //printf("%s", folder);
    sprintf(snapname, "%ssnapshot.csv.%03d", folder, snapnumber);
    snapfile = fopen(snapname, "w");
    fprintf(snapfile,"X,Y,Z,VX,VY,VZ,ID,T\n");
    for (n=0; n<ngals; n++){
        fprintf(snapfile,"%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f\n",
                gal[n].pos[0],
                gal[n].pos[1],
                gal[n].pos[2],
                gal[n].vel[0],
                gal[n].vel[1],
                gal[n].vel[2],
                n,
                t);
    }
    fclose(snapfile);
}
