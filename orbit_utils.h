/* ---------- coordinate conversion ---------- */
void convert(double *xtemp, double *vtemp, double *dsuntemp, double *vrsuntemp, double *vrtemp, double *ltemp, double *btemp, double *lcosbtemp, double *RAtemp, double *DECtemp, double *mu_alphatemp, double *mu_alphacosdeltatemp, double *mu_deltatemp, double *mutemp, double *PAtemp, int coordtype_coco, int vcoordtype, int radiococo, double vLSRtemp, double rgalsun){

    double x,y,z;        //galactic coordinates Input: [pc] Output: [kpc] (converted to kpc below)
    double xsun = -rgalsun/1000.0;//galactocentric distance of sun [kpc]
    double dsun_coco;           //heliocentric radial distance Input: [pc] Output: [kpc] (converted to kpc below)
    double dxy;             //heliocentric distance in xy-plane [kpc]
    double vx,vy,vz;        //cluster velocity [km/s]
    double vrx,vry,vrz;     //cluster radial velocity in 3d coordinates [km/s]
    double vtx,vty,vtz;     //cluster tangential velocity in 3d coordinates [km/s]
    double T[3][3], A[3][3], B[3][3];
    double TI[3][3];
    double RArad, DECrad;
    double brad, lrad, lcosbrad;
    double mu_alphacosdelta_coco, mu_delta_coco, mu, PArad;
    double vr_coco, vrsun;
    double RAENP = 0.0, DECENP = PI/2.0;  //equatorial coordinates of equatorial north pole
    double xENP, yENP, zENP, dxyENP; //cartesian vector pointing to the equatorial north pole
    double bENP, lENP; //galactic coordinates of ENP
    double FAK;
    double xdelta, ydelta, zdelta;
    double nx, ny, nz, dvt;
    double vrLSR, vrGSR;


    //transformation matrix equ. -> gal. from Johnson & Soderblom (1987)
/*  double t,d,a;
    double detT;
    t = PAGNP/360.0*2.0*PI;
    d = deltaGNP/360.0*2.0*PI;
    a = alphaGNP/360.0*2.0*PI;

    T[0][0] = -cos(t)*sin(d)*cos(a)-sin(t)*sin(a);
    T[0][1] = -cos(t)*sin(d)*sin(a)+sin(t)*cos(a);
    T[0][2] = cos(t)*cos(d);

    T[1][0] = -sin(t)*sin(d)*cos(a)+cos(t)*sin(a);
    T[1][1] = -sin(t)*sin(d)*sin(a)-cos(t)*cos(a);
    T[1][2] = sin(t)*cos(d);

    T[2][0] = cos(d)*cos(a);
    T[2][1] = cos(d)*sin(a);
    T[2][2] = sin(d);

    //invert matrix T in the most general way
    detT = T[0][0]*T[1][1]*T[2][2] + T[1][0]*T[2][1]*T[0][2] + T[2][0]*T[0][1]*T[1][2] - T[0][0]*T[1][2]*T[2][1] - T[1][0]*T[2][2]*T[0][1] - T[2][0]*T[0][2]*T[1][1];

    TI[0][0] = (T[1][1]*T[2][2]-T[1][2]*T[2][1])/detT;
    TI[1][0] = (T[1][2]*T[2][0]-T[1][0]*T[2][2])/detT;
    TI[2][0] = (T[1][0]*T[2][1]-T[1][1]*T[2][0])/detT;

    TI[0][1] = (T[0][2]*T[2][1]-T[0][1]*T[2][2])/detT;
    TI[1][1] = (T[0][0]*T[2][2]-T[2][0]*T[0][2])/detT;
    TI[2][1] = (T[0][1]*T[2][0]-T[0][0]*T[2][1])/detT;

    TI[0][2] = (T[0][1]*T[1][2]-T[0][2]*T[1][1])/detT;
    TI[1][2] = (T[0][2]*T[1][0]-T[0][0]*T[1][2])/detT;
    TI[2][2] = (T[0][0]*T[1][1]-T[0][1]*T[1][0])/detT;
*/

    //use result right away, careful when changing epochs though
    T[0][0] = -0.0548765333890783;
    T[0][1] = -0.8734366584039325;
    T[0][2] = -0.4838356847519307;
    T[1][0] = 0.4941106597543148;
    T[1][1] = -0.4448303583524283;
    T[1][2] = 0.7469809958795511;
    T[2][0] = -0.8676653859641706;
    T[2][1] = -0.1980766418440660;
    T[2][2] = 0.4559851115501738;

    TI[0][0] = -0.0548765333890783;
    TI[0][1] = 0.4941106597543147;
    TI[0][2] = -0.8676653859641704;
    TI[1][0] = -0.8734366584039324;
    TI[1][1] = -0.4448303583524283;
    TI[1][2] = -0.1980766418440659;
    TI[2][0] = -0.4838356847519305;
    TI[2][1] = 0.7469809958795510;
    TI[2][2] = 0.4559851115501737;


    //convert to kpc
    x = xtemp[0]/1000.0;
    y = xtemp[1]/1000.0;
    z = xtemp[2]/1000.0;

    dsun_coco = *dsuntemp/1000.0;

    vx = vtemp[0];
    vy = vtemp[1];
    vz = vtemp[2];

    vr_coco = *vrtemp;
    vrsun = *vrsuntemp;

    //convert to radians
    DECrad = *DECtemp/360.0*2.0*PI;
    RArad = *RAtemp/360.0*2.0*PI;
    PArad = *PAtemp/360.0*2.0*PI;


    //get the galactic coordinates first
    if (coordtype_coco == 1) {
        if (radiococo) printf("\nConverting equatorial to galactic coordinates using the transformation matrix:\n");
        if (radiococo) printf("%f\t%f\t%f\n",T[0][0],T[0][1],T[0][2]);
        if (radiococo) printf("%f\t%f\t%f\n",T[1][0],T[1][1],T[1][2]);
        if (radiococo) printf("%f\t%f\t%f\n",T[2][0],T[2][1],T[2][2]);

        brad = asin(T[2][0]*cos(DECrad)*cos(RArad) + T[2][1]*cos(DECrad)*sin(RArad) + T[2][2]*sin(DECrad));
        if (asin((T[1][0]*cos(DECrad)*cos(RArad) + T[1][1]*cos(DECrad)*sin(RArad) + T[1][2]*sin(DECrad))/cos(brad))>=0.0) {
            lrad = acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad));
        } else {
            lrad = 2.0*PI-acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad));
        }
        lcosbrad = lrad*cos(brad);
    } else if (coordtype_coco == 2) {
        brad = *btemp/360.0*2.0*PI;
        if (*ltemp) {
            lrad = *ltemp/360.0*2.0*PI;
            lcosbrad = lrad*cos(brad);
        } else if (*lcosbtemp) {
            lcosbrad = *lcosbtemp/360.0*2.0*PI;
            lrad = lcosbrad/cos(brad);
        }
    } else if (coordtype_coco == 3) {
        if (y >= 0.0)
            lrad = acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
        else
            lrad = 2.0*Pi-acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
        brad =  atan(z/sqrt(pow(x-xsun,2)+y*y));
        lcosbrad = lrad*cos(brad);
    }


    //get 3d position of cluster [kpc] from galactic coordinates
    if (coordtype_coco < 3) {
        z = sin(brad)*dsun_coco;
        dxy = sqrt(dsun_coco*dsun_coco-z*z);
        x = cos(lrad)*dxy + xsun;
        y = sin(lrad)*dxy;
    } else if (coordtype_coco == 3) {
        dsun_coco = sqrt(pow(x-xsun,2)+y*y+z*z);
    }


    //finally get the equatorial coordinates from galactic coordinates
    if (coordtype_coco > 1) {
        if (radiococo) printf("\nConverting galactic to equatorial coordinates using the transformation matrix:\n");
        if (radiococo) printf("%f\t%f\t%f\n",TI[0][0],TI[0][1],TI[0][2]);
        if (radiococo) printf("%f\t%f\t%f\n",TI[1][0],TI[1][1],TI[1][2]);
        if (radiococo) printf("%f\t%f\t%f\n",TI[2][0],TI[2][1],TI[2][2]);

        if (radiococo) {
            //unit matrix B = T * TI
            B[0][0] = T[0][0]*TI[0][0] + T[0][1]*TI[1][0] + T[0][2]*TI[2][0];
            B[0][1] = T[0][0]*TI[0][1] + T[0][1]*TI[1][1] + T[0][2]*TI[2][1];
            B[0][2] = T[0][0]*TI[0][2] + T[0][1]*TI[1][2] + T[0][2]*TI[2][2];

            B[1][0] = T[1][0]*TI[0][0] + T[1][1]*TI[1][0] + T[1][2]*TI[2][0];
            B[1][1] = T[1][0]*TI[0][1] + T[1][1]*TI[1][1] + T[1][2]*TI[2][1];
            B[1][2] = T[1][0]*TI[0][2] + T[1][1]*TI[1][2] + T[1][2]*TI[2][2];

            B[2][0] = T[2][0]*TI[0][0] + T[2][1]*TI[1][0] + T[2][2]*TI[2][0];
            B[2][1] = T[2][0]*TI[0][1] + T[2][1]*TI[1][1] + T[2][2]*TI[2][1];
            B[2][2] = T[2][0]*TI[0][2] + T[2][1]*TI[1][2] + T[2][2]*TI[2][2];

            printf("\nCalculating T*T^{-1} = 1 for consistency check:\n");
            printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
            printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
            printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
        }

        DECrad = asin(TI[2][0]*cos(brad)*cos(lrad)+TI[2][1]*cos(brad)*sin(lrad)+TI[2][2]*sin(brad));
        if (asin((TI[1][0]*cos(brad)*cos(lrad) + TI[1][1]*cos(brad)*sin(lrad) + TI[1][2]*sin(brad))/cos(DECrad))>=0.0) {
            RArad = acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad));
        } else {
            RArad = 2.0*PI-acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad));
        }
    }



    //get tangential velocity in [km/s] from different sets of velocity-measurement types

    //get coordinates of equatorial north pole on great circle
    bENP = asin(T[2][0]*cos(DECENP)*cos(RAENP) + T[2][1]*cos(DECENP)*sin(RAENP) + T[2][2]*sin(DECENP));
    if (asin((T[1][0]*cos(DECENP)*cos(RAENP) + T[1][1]*cos(DECENP)*sin(RAENP) + T[1][2]*sin(DECENP))/cos(bENP))>=0.0) {
        lENP = acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP));
    } else {
        lENP = 2.0*PI-acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP));
    }
    if (radiococo) printf("\nCoordinates of equatorial north pole:\n");
    if (radiococo) printf("bENP = %f\tlENP = %f\n", bENP, lENP);
    zENP = sin(bENP)*dsun_coco;
    dxyENP = sqrt(dsun_coco*dsun_coco-zENP*zENP);
    xENP = cos(lENP)*dxyENP + xsun;
    yENP = sin(lENP)*dxyENP;
    if (radiococo) printf("xENP = %f\tyENP = %f\tzENP = %f\n", xENP, yENP, zENP);


    if (vcoordtype == 1) {

        //get radial velocity in 3d coordinates [km/s]
        vrx = (x - xsun)/dsun_coco*vrsun;
        vry = y/dsun_coco*vrsun;
        vrz = z/dsun_coco*vrsun;
        if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
        if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));

        //convert to km/s
        mu = *mutemp*dsun_coco*4.74057;

        //compute proper motion components
        mu_alphacosdelta_coco = mu*sin(PArad);
        mu_delta_coco = mu*cos(PArad);

        A[0][0] = cos(RArad)*cos(DECrad);
        A[0][1] = -sin(RArad);
        A[0][2] = -cos(RArad)*sin(DECrad);

        A[1][0] = sin(RArad)*cos(DECrad);
        A[1][1] = cos(RArad);
        A[1][2] = -sin(RArad)*sin(DECrad);

        A[2][0] = sin(DECrad);
        A[2][1] = 0.0;
        A[2][2] = cos(DECrad);

        //printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
        //printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
        //printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);

        //B = T * A
        B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
        B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
        B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];

        B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
        B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
        B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];

        B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
        B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
        B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];

        //printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
        //printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
        //printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);

        vx = vrsun*B[0][0] + mu_alphacosdelta_coco*B[0][1] + mu_delta_coco*B[0][2] +vxsun;
        vy = vrsun*B[1][0] + mu_alphacosdelta_coco*B[1][1] + mu_delta_coco*B[1][2] +vysun+vLSRtemp;
        vz = vrsun*B[2][0] + mu_alphacosdelta_coco*B[2][1] + mu_delta_coco*B[2][2] +vzsun;


        if (radiococo) printf("\nCartesian velocity:\n");
        if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
        if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));

    } else if (vcoordtype == 2) {

        //get radial velocity in 3d coordinates [km/s]
        vrx = (x - xsun)/dsun_coco*vrsun;
        vry = y/dsun_coco*vrsun;
        vrz = z/dsun_coco*vrsun;
        if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
        if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));

        if (*mu_alphatemp) *mu_alphacosdeltatemp = *mu_alphatemp*cos(DECrad);
        else if (*mu_alphacosdeltatemp) *mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);

        //convert to km/s
        mu_alphacosdelta_coco = *mu_alphacosdeltatemp*dsun_coco*4.74057;
        mu_delta_coco = *mu_deltatemp*dsun_coco*4.74057;
        mu = sqrt(mu_alphacosdelta_coco*mu_alphacosdelta_coco+mu_delta_coco*mu_delta_coco);

        A[0][0] = cos(RArad)*cos(DECrad);
        A[0][1] = -sin(RArad);
        A[0][2] = -cos(RArad)*sin(DECrad);

        A[1][0] = sin(RArad)*cos(DECrad);
        A[1][1] = cos(RArad);
        A[1][2] = -sin(RArad)*sin(DECrad);

        A[2][0] = sin(DECrad);
        A[2][1] = 0.0;
        A[2][2] = cos(DECrad);

        //printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
        //printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
        //printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);

        //B = T * A
        B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
        B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
        B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];

        B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
        B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
        B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];

        B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
        B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
        B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];

        //printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
        //printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
        //printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);

        vx = vrsun*B[0][0] + mu_alphacosdelta_coco*B[0][1] + mu_delta_coco*B[0][2] +vxsun;
        vy = vrsun*B[1][0] + mu_alphacosdelta_coco*B[1][1] + mu_delta_coco*B[1][2] +vysun+vLSRtemp;
        vz = vrsun*B[2][0] + mu_alphacosdelta_coco*B[2][1] + mu_delta_coco*B[2][2] +vzsun;

        if (radiococo) printf("\nCartesian velocity:\n");
        if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
        if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));

        //get position angle of proper motion

        //heliocentric transverse velocity
        vtx = vx-vxsun-vrx;
        vty = vy-vysun-vLSRtemp-vry;
        vtz = vz-vzsun-vrz;
        if (radiococo) printf("\nTransverse velocity:\n");
        if (radiococo) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));

        //get tangential vector pointing to ENP
        FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
        xdelta = FAK*(x-xsun)+(xENP-xsun);
        ydelta = FAK*y+yENP;
        zdelta = FAK*z+zENP;

        //determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
        nx = y*zENP-z*yENP;
        ny = z*(xENP-xsun)-(x-xsun)*zENP;
        nz = (x-xsun)*yENP-y*(xENP-xsun);
        dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;

        //get position angle of proper motion with respect to tangential vector pointing to ENP
        if (dvt <= 0)
            PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
        else
            PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));

        if (radiococo) printf("\nProper motion and position angle:\n");
        if (radiococo) printf("mu = %f\tPA = %f\n", mu, PArad);

    } else if (vcoordtype == 3) {

        if (radiococo) printf("\nCartesian velocity:\n");
        if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
        if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));

        //heliocentric radial velocity
        vrsun = ((vx-vxsun)*(x-xsun)+(vy-vysun-vLSRtemp)*y+(vz-vzsun)*z)/sqrt(pow(x-xsun,2)+y*y+z*z);

        //get radial velocity in 3d coordinates [km/s]
        vrx = (x - xsun)/dsun_coco*vrsun;
        vry = y/dsun_coco*vrsun;
        vrz = z/dsun_coco*vrsun;
        if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
        if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));

        //get position angle of proper motion

        //heliocentric transverse velocity
        vtx = vx-vxsun-vrx;
        vty = vy-vysun-vLSRtemp-vry;
        vtz = vz-vzsun-vrz;
        if (radiococo) printf("\nTransverse velocity:\n");
        if (radiococo) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));

        //get tangential vector pointing to ENP
        FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
        xdelta = FAK*(x-xsun)+(xENP-xsun);
        ydelta = FAK*y+yENP;
        zdelta = FAK*z+zENP;

        //determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
        nx = y*zENP-z*yENP;
        ny = z*(xENP-xsun)-(x-xsun)*zENP;
        nz = (x-xsun)*yENP-y*(xENP-xsun);
        dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;

        //get position angle of proper motion with respect to tangential vector pointing to ENP
        if (dvt <= 0)
            PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
        else
            PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));

        if (radiococo) printf("\nProper motion and position angle:\n");
        if (radiococo) printf("mu = %f\tPA = %f\n", mu, PArad);

        mu = sqrt(vtx*vtx+vty*vty+vtz*vtz);
        mu_delta_coco = mu*cos(PArad);
        mu_alphacosdelta_coco = mu*sin(PArad);

    }


    if (radiococo) printf("\nProper motion:\n");
    if (radiococo) printf("mu_alphacosdelta  = %f\tmu_delta = %f\tmu = %f [km/s]\t PA = %f\n", mu_alphacosdelta_coco, mu_delta_coco, mu, PArad);

    vr_coco = (vx*(x-xsun)+vy*y+vz*z)/sqrt(pow(x-xsun,2)+y*y+z*z);
    if (radiococo) printf("\nRadial velocity:\n");
    if (radiococo) printf("vr = %.3f\tvr = %.3f (heliocentric) [km/s]\n", vr_coco, vrsun);

    //consistency check with formula for GSR radial velocity from script of Steven Majewski
    if (radiococo) vrLSR = vrsun + (vxsun*cos(brad)*cos(lrad)+vysun*cos(brad)*sin(lrad)+vzsun*sin(brad));
    if (radiococo) vrGSR = vrLSR + vLSRtemp*cos(brad)*sin(lrad);
    if (radiococo) printf("\nConsistency check with formula for Galactic standard of rest (GSR) radial velocity (should be equal to vr):\n");
    if (radiococo) printf("vr_LSR = %f\tvr_GSR = %.3f [km/s]\n", vrLSR, vrGSR);



    //convert back to input units and write to output
    *xtemp = 1000.0*x;
    *(xtemp+1) = 1000.0*y;
    *(xtemp+2) = 1000.0*z;

    *dsuntemp = 1000.0*dsun_coco;

    *vtemp = vx;
    *(vtemp+1) = vy;
    *(vtemp+2) = vz;

    *vrsuntemp = vrsun;
    *vrtemp = vr_coco;

    *DECtemp = DECrad*180.0/PI;
    *RAtemp = RArad*180.0/PI;

    *btemp = brad*180.0/PI;
    *ltemp = lrad*180.0/PI;
    *lcosbtemp = *ltemp*cos(brad);

    *mutemp = mu/(dsun_coco*4.74057);
    *PAtemp = PArad*180.0/PI;
    *mu_deltatemp = mu_delta_coco/(dsun_coco*4.74057);
    *mu_alphacosdeltatemp = mu_alphacosdelta_coco/(dsun_coco*4.74057);
    *mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);

}

/* ---------- gaussian distribution ---------- */
double get_gauss(void){
    double random[2],p,q;
    do {
        random[0] = 2.0*drand48()-1.0;
        random[1] = 2.0*drand48()-1.0;
        q = random[0]*random[0]+random[1]*random[1];
    } while (q>1.0);

    p = sqrt(-2.0*log(q)/q);
    return random[0]*p;

}

/* ---------- sorting functions ------------ */
void shellsort(double **array, int N, int k) {//largest up
    int i,j,n,o;
    N = N-1;
    double swap[k];
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            for (o=0; o<k; o++) swap[o] = array[i][o];
            for (j = i; ((j >= n) && (array[j-n][0] < swap[0])); j -= n) {
                for (o=0; o<k; o++) array[j][o] = array[j-n][o];
            }
            for (o=0; o<k; o++) array[j][o] = swap[o];
        }
    }
}

void shellsort_reverse(double **array, int N, int k) {//smallest up
    int i,j,o,n;
    N = N-1;
    double swap[k];
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            for (o=0; o<k; o++) swap[o] = array[i][o];
            for (j = i; ((j >= n) && (array[j-n][0] > swap[0])); j -= n) {
                for (o=0; o<k; o++) array[j][o] = array[j-n][o];
            }
            for (o=0; o<k; o++) array[j][o] = swap[o];
        }
    }
}

void shellsort_1d(double *array, int N) {//largest up
    int i,j,n;
    N = N-1;
    double swap;
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            swap = array[i];
            for (j = i; ((j >= n) && (array[j-n] < swap)); j -= n) {
                array[j] = array[j-n];
            }
            array[j] = swap;
        }
    }
}

void shellsort_reverse_1d(double *array, int N) {//smallest up
    int i,j,n;
    N = N-1;
    double swap;
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            swap = array[i];
            for (j = i; ((j >= n) && (array[j-n] > swap)); j -= n) {
                array[j] = array[j-n];
            }
            array[j] = swap;
        }
    }
}

