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
    // Note: no input error checking done here. Do that in python wrapper.    

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
    double tfuture = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "tfuture"));
    double dt0 = PyFloat_AsDouble(PyDict_GetItemString(input_parameters, "dt0"));
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
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal, parameters, sign);

    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
        sign = 1.0;
        dtoutt = dtout;
        tmax = tfuture;
        err = rk4_drv(&t, tmax, dtoutt, dt0, mdiff, gal, parameters, sign);
        printf("%d", err); 
    }
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
	int dir = -1;
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
		    *t+=dt;
		    //update the test particles here
                    for (n=0; n<ngals; n++){
		        for (k=0;k<3;k++) {
			    (*(gal+n)).pos[k]=(*(gal+n)).post[k];/* If yes -> continue and double step size */
			    (*(gal+n)).vel[k]=(*(gal+n)).velt[k];
		        }
		    }
		    if (VARIABLE_TIMESTEPS) {
                        dt = dt*2.0;
                    }
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
                if (snapshot){
                    write_snapshot(parameters, gal, *t, snapnum);
                    snapnum += 1;
                }
                if (orbit_stats) {
                    //record orbit stats here (pericenters and apocenters etc.)
                    // could also do it above but I think that would be add too many operations per timestep
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


/* ----------- force ----------- */
void getforce(double *x, double *v, double *a, struct Params parameters, struct Gal gal){
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

        if (parameters.Mhalo > 0.0){

            r3 = sqrt(pow(*x, 2) + pow(*(x+1), 2) + pow(*(x+2)/q_halo, 2));

            // NFW
            if (NFW) {
                ax += -2.0*G*parameters.Mhalo * (log(1.0 + r3/parameters.r_halo)/r3 -
                                                1.0/(parameters.r_halo+r3)) * *x/pow(r3, 2);
                ay += -2.0*G*parameters.Mhalo * (log(1.0 + r3/parameters.r_halo)/r3 -
                                                1.0/(parameters.r_halo+r3)) * *(x+1)/pow(r3, 2);
                az += -2.0*G*parameters.Mhalo * (log(1.0 + r3/parameters.r_halo)/r3 -
                                                1.0/(parameters.r_halo+r3)) * *(x+2)/
                                                pow(parameters.q_halo*parameters.q_halo*r3, 2);

            } else {
            //Dehnen Halo
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
                // Coulomb logarithm using Dehnen mass
              //  double coulomb = r3/(1.4*3000); 

                double coulomb = r3*vr*vr/(G*gal.mhalo*
                                pow(r3/(r3+gal.r_halo), 3-gal.gamma));  
                // This assumes that vcirc == sqrt(3)*vdisp and that the dispersion is measured at the scale radius
                double X;
                if (NFW) {
                    X = 1.0/(sqrt(2.0*log(2.0)-1))*sqrt(parameters.r_halo/(parameters.Mhalo*G))*vr;
                } else {
                    X = 2*sqrt(parameters.r_halo/(parameters.Mhalo*G))*vr;
                }
                //double sigma = sqrt(G*parameters.Mhalo*r3*pow(parameters.r_halo + r3, 3)*
                //        (-(25.0*pow(parameters.r_halo, 3) + 52.0*pow(parameters.r_halo, 2)*r3 + 
                //            42.0*parameters.r_halo*pow(r3, 2) + 12.0*pow(r3, 3))/
                //            (12.0*pow(parameters.r_halo, 4)*pow(parameters.r_halo + r3, 4)) + 
                //            log((parameters.r_halo + r3)/r3)/pow(parameters.r_halo, 5)));
                //double X = vr/(sqrt(2)*sigma); 
     
                // This assumes that vcirc == sqrt(2)*vdisp and that the mass profile is spherical
                double density;
                if (NFW) {
                    // where rho0 = Mhalo/(2Pi*a^3)
                    density = parameters.Mhalo*parameters.r_halo/(2.0*Pi)/
                                        (r3*pow(r3+parameters.r_halo, 2));             
                } else {
                    density = (3 - parameters.gamma)*parameters.Mhalo/(4*Pi)*
                            parameters.r_halo/(pow(r3, parameters.gamma)*pow(r3 + parameters.r_halo,
                                                                            4-parameters.gamma));
                }

                ax += -4.0*Pi*G*G*gal.mhalo*density*log(coulomb)*
                            (1.0/6.0)*(erf(X) - 2*X/sqrt(Pi)*exp(-X*X))* *v/pow(vr, 3);       
                ay += -4.0*Pi*G*G*gal.mhalo*density*log(coulomb)*
                            (1.0/6.0)*(erf(X) - 2*X/sqrt(Pi)*exp(-X*X))* *(v+1)/pow(vr, 3);       
                az += -4.0*Pi*G*G*gal.mhalo*density*log(coulomb)*
                            (1.0/6.0)*(erf(X) - 2*X/sqrt(Pi)*exp(-X*X))* *(v+2)/pow(vr, 3);       
            }
        }
	*(a+0) = ax;
	*(a+1) = ay;
	*(a+2) = az;
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

    double coulomb = 0.0;
    double X = 0.0;
    double density = 0.0; 
    double softening = 1.4*100;  //epsilon = 0.1kpc
    double sigma = 0.0;

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


            //Dehnen Halo
            r = sqrt(pow(*x - gal[i].pos[0], 2) +
                     pow(*(x+1) - gal[i].pos[1], 2) +
                     pow(*(x+2) - gal[i].pos[2], 2));
            ax += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*x - gal[i].pos[0]);
            ay += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*(x+1) - gal[i].pos[1]);
            
            az += -G*gal[i].mhalo * pow(r/(gal[i].r_halo + r), -gal[i].gamma)/
                                    pow(gal[i].r_halo + r, 3) * (*(x+2) - gal[i].pos[2]);

            // dynamical friction
            if (DYNAMICALFRICTION_COMPANION && r > softening) {

                //relative velocity
                vx = (*v - gal[i].vel[0]);
                vy = (*(v+1) - gal[i].vel[1]);
                vz = (*(v+2) - gal[i].vel[2]);
                vr = sqrt(vx*vx + vy*vy + vz*vz); 
                // Coulomb logarithm using Dehnen mass
                // use the r/1.6*epsilon formulation
               // coulomb = r/(softening); // epsilon = 0.1 kpc

                coulomb = r*vr*vr/(G*gal[gal_num].mhalo*
                               pow(r/(r+gal[gal_num].r_halo), 3-gal[gal_num].gamma));  
                // This assumes that vcirc == sqrt(2)*vdisp and that the mass profile is spherical
            //  X = 2*sqrt(gal[i].r_halo/(gal[i].mhalo*G))*vr;
                // Hernquist 1D velocity dispersion - CForm from Mathematica
                sigma = 3.0*sqrt(G*gal[i].mhalo*r*pow(gal[i].r_halo + r, 3)*
                        (-(25.0*pow(gal[i].r_halo, 3) + 52.0*pow(gal[i].r_halo, 2)*r + 
                            42.0*gal[i].r_halo*pow(r, 2) + 12.0*pow(r, 3))/
                            (12.0*pow(gal[i].r_halo, 4)*pow(gal[i].r_halo + r, 4)) + 
                            log((gal[i].r_halo + r)/r)/pow(gal[i].r_halo, 5)));
                X = vr/(sqrt(2)*sigma); 
                
                density = (3 - gal[i].gamma)*gal[i].mhalo/(4*Pi)*
                            gal[i].r_halo/(pow(r, gal[i].gamma)*pow(r + gal[i].r_halo, 4-gal[i].gamma));
    
                ax += -4.0*Pi*G*G*gal[gal_num].mhalo*density*log(coulomb)*
                        (1.0/1.0)*(erf(X) - 2*X/sqrt(Pi)*exp(-X*X))*vx/pow(vr, 3);       

                ay += -4.0*Pi*G*G*gal[gal_num].mhalo*density*log(coulomb)*
                        (1.0/1.0)*(erf(X) - 2*X/sqrt(Pi)*exp(-X*X))*vy/pow(vr, 3);       

                az += -4.0*Pi*G*G*gal[gal_num].mhalo*density*log(coulomb)*
                        (1.0/1.0)*(erf(X) - 2*X/sqrt(Pi)*exp(-X*X))*vz/pow(vr, 3);       

            }

        }
    }
    *(a+0) += ax;
    *(a+1) += ay;
    *(a+2) += az;

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
    else {
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


void write_snapshot(struct Params parameters, struct Gal *gal, double t, int snapnumber){
    int n;
    int ngals = parameters.ngals;
    char *folder = parameters.outputdir;
    FILE *snapfile;
    char snapname[50];
    double acc0[3];
    //printf("%s", folder);
    sprintf(snapname, "%ssnapshot.csv.%03d", folder, snapnumber);
    snapfile = fopen(snapname, "w");
    fprintf(snapfile,"X,Y,Z,VX,VY,VZ,AX,AY,AZ,T\n");
    for (n=0; n<ngals; n++){
        getforce_gals(gal[n].pos, gal[n].vel, acc0, n, gal, parameters);
        fprintf(snapfile,"%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f\n",
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
