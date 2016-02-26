#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include <sys/time.h>
#include <unistd.h>
#include <omp.h>
#include <getparam.h> /* make library in directory libminimo first */
#include "muesli.h"



int main (int argc, char **argv) {
    //***Allocate***
    t0 = malloc(BMAX*sizeof(double));
    mass = malloc(BMAX*sizeof(double));
    x = malloc(BMAX*sizeof(double));
    y  = malloc(BMAX*sizeof(double));
    z = malloc(BMAX*sizeof(double));
    vx = malloc(BMAX*sizeof(double));
    vy = malloc(BMAX*sizeof(double));
    vz = malloc(BMAX*sizeof(double));
    phif = malloc(BMAX*sizeof(double));
    mgc = malloc(BMAX*sizeof(double));
    rh0 = malloc(BMAX*sizeof(double));
    rh = malloc(BMAX*sizeof(double));
    laddx = malloc(BMAX*sizeof(double));
    laddy = malloc(BMAX*sizeof(double));
    laddz = malloc(BMAX*sizeof(double));
    tboxo = malloc(BMAX*sizeof(double));
    tboxn = malloc(BMAX*sizeof(double));
    ts = malloc(BMAX*sizeof(double));
    box = malloc(BMAX*sizeof(int));
    boy = malloc(BMAX*sizeof(int));
    boz = malloc(BMAX*sizeof(int));
    boxorb = malloc(BMAX*sizeof(int));
    grid_cell_nr = malloc(BMAX*sizeof(int));
    va =  malloc(BMAX*sizeof(double));            //GC Mass loss by Relaxation
    vp =  malloc(BMAX*sizeof(double));
    mgct0 = malloc(BMAX*sizeof(double));
    tgc = malloc(BMAX*sizeof(double));
    rperi =  malloc(BMAX*sizeof(double));
    rapo =  malloc(BMAX*sizeof(double));
    tdisso =  malloc(BMAX*sizeof(double));
    
    
    
    
    get_host(argv[0]);
    get_coeff0();
    initparam(argv, defv);
    
    if (streq(getparam("restart"),"Y") || streq(getparam("restart"),"y")) {
        
        read_restart();
        printf("\nRestarting old run from T = %lf\n",tout);
        
    } else {
        
        tdis      = getdparam("tdis");
        gcevolu   = getdparam("gcevolu");
        gasexp    = getdparam("gasexp");
        scalemass = getdparam("scalemass");
        scalesize = getdparam("scalesize");
        xcrit     = getdparam("xcrit");
        dfmethod  = getdparam("dfmethod");
        tend      = getdparam("tend");
        rcap      = getdparam("rcap");
        msmbh     = getdparam("msmbh");
        msmbh0    = getdparam("msmbh");
        read_nbody();
        
        printf("\n---------------------- Start integration -----------------------\n");
        
        tout   = 0.0;
        twrite = 0.0;
        
    }
    
    
    
    printf("\nNumber of stars:     %i\n",nbody);
    printf("End of simulation:   %lf\n",tend);
    
    if (tdis == 1) {
        printf("Van Hoerner Method for TD %i\n",tdis);}
    else if (tdis == 2) {
        printf("Kuepper Method for TD %i\n",tdis);}
    else if (tdis == 0) {
        printf("No GC Disruption %i\n",tdis);}
    
    if (gasexp == 1) {
        printf("Gas expulsion on %i\n",gasexp);}
    else {
        printf("Gas expulsion off %i\n",gasexp);}
    
    
    if (gcevolu == 1) {
        printf("GC mass evolution on %i\n",gcevolu);}
    else {
        printf("GC mass evolution off %i\n",gcevolu);}
    
    if (dfmethod == 1) {
        printf("Use standard Chandrasekhar formula for DF %i\n",dfmethod);}
    else if (dfmethod == 2) {
        printf("Use triaxial generalization of Chandrasekhar formular for DF %i\n",dfmethod);}
    else {
        printf("Dynamical Friction off %i\n",dfmethod);}
    
    printf("Scalefactor mass [M_sun]:     %le\n",scalemass);
    printf("Scalefactor size [pc]:     %le\n",scalesize);
    printf("BH capture radius:  %le\n",rcap);
    printf("SMBH mass fraction:  %lf\n\n",msmbh);
    
    mt_drv();
    
    printf("\n----------------------- End integration ------------------------\n");
    
    if (rcap>0.0) printf("\nTotal number of tidal disruptions: %i\n",ncap);
    
    
    //***Deallocate***
    free(t0);
    free(mass);
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(phif);
    free(mgc);
    free(rh0);
    free(rh);
    free(laddx);
    free(laddy);
    free(laddz);
    free(tboxo);
    free(tboxn);
    free(ts);
    free(box);
    free(boy);
    free(boz);
    free(boxorb);
    free(grid_cell_nr);
    free(va);                   //GC Mass loss by Relaxation
    free(vp);
    free(mgct0);
    free(tgc);
    free(rperi);
    free(rapo);
    free(tdisso);
    return 0;
    
}


void get_host(char *pname) {
    
    int i,size=100;
    char *name,*dirn;
    time_t timep;
    
    printf("\n\n----------------------------------------------------------------\n");
    
    name = (char *) malloc (size*sizeof(char));
    dirn = (char *) malloc (size*sizeof(char));
    getcwd(dirn,size);
    printf("ISTAMP - This is program: %s/%s\n",dirn,pname+2);
    i = gethostname(name, size);
    
#pragma omp parallel
#pragma omp master
    {
        printf("ISTAMP - running on %s with %d CPU core(s)\n",name,omp_get_num_threads( ));
    }
    
    printf("ISTAMP - PID: %i\n",getpid());
    time(&timep);
    printf("ISTAMP - Program (re-)started %s",ctime(&timep));
    printf("----------------------------------------------------------------\n");
    
}


void read_nbody() {
    
    FILE *instr;
    char *fname=calloc (100,sizeof(char));
    char *dname=calloc (100,sizeof(char));
    int i;
    
    strcpy (dname,getparam("in"));
    strcpy (fname,dname);
    strcat (fname,".NBODY");
    
    printf("\nReading Input-File: %s\n",fname);
    
    if ((instr = fopen(fname,"r"))==NULL) {
        printf("Input File not found !\n");
        exit(-1);
    }
    
    while (fscanf(instr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",(mass+nbody),(mgc+nbody),(rh0+nbody),(x+nbody),(y+nbody),(z+nbody),(vx+nbody),(vy+nbody),(vz+nbody))!=EOF) nbody++;
    
    fclose(instr);
    
    printf("\nFinished reading %i stars !\n\n",nbody);
    
    if (scalemass == 0.0 || scalesize == 0.0) {
        exit(-1);
    }
    
    
    for (i=0;i<nbody;i++)  t0[i] = 0.0;
    
    strcpy(fname,dname);
    strcat (fname,".POS");
    instr = fopen(fname,"w");
    
    fclose(instr);
    
}


void write_nbody(double time) {
    
    FILE *dat;
    char *fname=calloc (100,sizeof(char));
    strcpy (fname,getparam("in"));
    strcat (fname,".POS");
    dat=fopen(fname,"a");
    
    
    fwrite (&time,sizeof(double),1,dat);
    fwrite (&nbody,sizeof(int),1,dat);
    fwrite (mass,sizeof(double),nbody,dat);
    fwrite (mgc,sizeof(double),nbody,dat);
    fwrite (rh ,sizeof(double),nbody,dat);
    fwrite (x,sizeof(double),nbody,dat);
    fwrite (y,sizeof(double),nbody,dat);
    fwrite (z,sizeof(double),nbody,dat);
    fwrite (vx,sizeof(double),nbody,dat);
    fwrite (vy,sizeof(double),nbody,dat);
    fwrite (vz,sizeof(double),nbody,dat);
    fwrite (phif,sizeof(double),nbody,dat);
    
    fclose(dat);
    
}


void write_restart() {
    
    FILE *dat;
    char *fname=calloc (100,sizeof(char));
    int i;
    
    // Before writing make a consistency check of the data
    
    for (i=0;i<nbody;i++) {
        if (isnan(mass[i])) { printf("Mass of body %i is NaN !\n",i); exit(-1); }
        if (isnan(mgc[i]))  { printf("Mass of GC %i is NaN !\n",i); exit(-1); }
        if (isnan(rh[i]))   { printf("Size of GC %i is NaN !\n",i); exit(-1); }
        if (isnan(x[i]))    { printf("x coordinate of body %i is NaN !\n",i); exit(-1); }
        if (isnan(y[i]))    { printf("y coordinate of body %i is NaN !\n",i); exit(-1); }
        if (isnan(z[i]))    { printf("z coordinate of body %i is NaN !\n",i); exit(-1); }
        if (isnan(vx[i]))   { printf("x velocity of body %i is NaN !\n",i); exit(-1); }
        if (isnan(vy[i]))   { printf("y velocity of body %i is NaN !\n",i); exit(-1); }
        if (isnan(vz[i]))   { printf("z velocity of body %i is NaN !\n",i); exit(-1); }
        if (isnan(t0[i]))   { printf("t0 of body %i is NaN !\n",i); exit(-1); }
        if (isnan(phif[i])) { printf("phif of body %i is NaN !\n",i); exit(-1); }
    }
    
    strcpy (fname,getparam("in"));
    strcat (fname,".RSTART");
    dat=fopen(fname,"w");
    
    fwrite (&nbody,sizeof(int),1,dat);
    fwrite (mass,sizeof(double),nbody,dat);
    fwrite (mgc,sizeof(double),nbody,dat);
    fwrite (rh,sizeof(double),nbody,dat);
    fwrite (x,sizeof(double),nbody,dat);
    fwrite (y,sizeof(double),nbody,dat);
    fwrite (z,sizeof(double),nbody,dat);
    fwrite (vx,sizeof(double),nbody,dat);
    fwrite (vy,sizeof(double),nbody,dat);
    fwrite (vz,sizeof(double),nbody,dat);
    fwrite (t0,sizeof(double),nbody,dat);
    fwrite (&rcap,sizeof(double),1,dat);
    fwrite (&tout,sizeof(double),1,dat);
    fwrite (&twrite,sizeof(double),1,dat);
    fwrite (&tend,sizeof(double),1,dat);
    fwrite (&msmbh,sizeof(double),1,dat);
    fwrite (&scalemass,sizeof(double),1,dat);
    fwrite (&scalesize,sizeof(double),1,dat);
    fwrite (&tdis,sizeof(int),1,dat);
    fwrite (&gcevolu,sizeof(int),1,dat);
    fwrite (&dfmethod,sizeof(int),1,dat);
    fclose(dat);
    
}


void read_restart() {
    
    FILE *dat;
    char *fname=calloc (100,sizeof(char));
    strcpy (fname,getparam("in"));
    strcat (fname,".RSTART");
    dat=fopen(fname,"r");
    
    fread (&nbody,sizeof(int),1,dat);
    fread (mass,sizeof(double),nbody,dat);
    fread (mgc,sizeof(double),nbody,dat);
    fread (rh,sizeof(double),nbody,dat);
    fread (x,sizeof(double),nbody,dat);
    fread (y,sizeof(double),nbody,dat);
    fread (z,sizeof(double),nbody,dat);
    fread (vx,sizeof(double),nbody,dat);
    fread (vy,sizeof(double),nbody,dat);
    fread (vz,sizeof(double),nbody,dat);
    fread (t0,sizeof(double),nbody,dat);
    fread (&rcap,sizeof(double),1,dat);
    fread (&tout,sizeof(double),1,dat);
    fread (&twrite,sizeof(double),1,dat);
    fread (&tend,sizeof(double),1,dat);
    fread (&msmbh,sizeof(double),1,dat);
    fwrite (&scalemass,sizeof(double),1,dat);
    fwrite (&scalesize,sizeof(double),1,dat);
    fwrite (&tdis,sizeof(int),1,dat);
    fwrite (&gcevolu,sizeof(int),1,dat);
    fwrite (&dfmethod,sizeof(int),1,dat);
    fclose(dat);
    
}


void get_coeff0() {
    
    int n,l,m;
    double fac1,dm0,kln;
    
    for (l=0;l<LMAX;l++) {
        for (n=0;n<NMAX;n++) {
            
            kln = (double) 0.5*n*(n+4.0*l+3.0)+(l+1.0)*(2.0*l+1.0);
            abar[l][n] = (double) -1.0/kln*pow(2.0,8.0*l+6.0)/(4.0*PI)*factrl(n)*(n+2.0*l+1.5)*exp(2.0*gammln(2.0*l+1.5)-gammln(1.0*n+4.0*l+3.0));
            
        }
    }
    
    fac1 = 1.0/(4.0*PI);
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            if (m==0) dm0=1.0;
            else dm0=0.0;
            
            narr[l][m]= fac1*(2.0*l+1.0)*(2.0-dm0)*factrl(l-m)/factrl(l+m);
        }
    }
    
}


void mt_drv() {
    
    int i,k,nstep=0, cancel_loop;
    double a[3],pt,nstepg;
    double t1,t2,tv=0.0,tn=0.0,ns0=0.0,vscale,eccen=0.0;
    double vecr[3],veca[3];
    double rvo,rvn,dis,dt,accf,fst;
    double lxo,lyo,lzo,lxn,lyn,lzn;
    double dilaton,chit,trelax0;
    double xttl[3],vttl[3];
    FILE* outputtest;
    outputtest=fopen("outputtest.txt","w");
    tscale = 14.94*pow(pow(scalesize,3)/scalemass,0.5);   //in Myr
    vscale = 0.06557*pow(scalemass/scalesize,0.5);        //in Kms^-1
    
    for (i=0;i<nbody;i++) {
        ts[i]       =0.0;
        tboxo[i]    = 0.0;
        tboxn[i]    = 0.0;
        box[i]      =   0;
        boy[i]      =   0;
        boz[i]      =   0;
        boxorb[i]   =   0;
        laddx[i]    = 0.0;
        laddy[i]    = 0.0;
        laddz[i]    = 0.0;
        rh[i]=rh0[i];
    }
    
    
    if (gcevolu == 0) tinit=0.0;
    else if (gcevolu == 1) tinit=0.0; //Clusters are not disrupted (and no DF) within first xxx Myr
    
    
    
    // Gas expulsion (Kroupa & Boily (2002))
    if (gasexp == 1) {
        for (i=0;i<nbody;i++) {
            if (mgc[i] > 1.0E-15) {
                fst = exp(-0.5*pow((log10(mgc[i])-4.0)/0.5,2.0));
                mgc[i] *= (0.5-0.4*fst);
            }
        }
    }
    
    
    
    com_correction();
    get_coeff(tout);
    get_pot(tout);
    //  Definitions req. for Relaxation in Tidal Fields
    for (i=0;i<nbody;i++) {
        vecr[0]     = x[i];
        vecr[1]     = y[i];
        vecr[2]     = z[i];
        get_forcettl(vecr, veca);
        va[i]       = sqrt(sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i])*sqrt(veca[0]*veca[0]+veca[1]*veca[1]+veca[2]*veca[2]));
        vp[i]       = va[i];
        mgct0[i]    = mgc[i];
        tgc[i]      = 0.0;
        rperi[i]    = 0.0;
        rapo[i]     = 0.0;
        if (mgc[i] > 100.0) {
            ns0 = mgct0[i]/0.64;
            tdisso[i] = 1.91*pow(ns0/(log(0.02*ns0)),0.75)*sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i])*scalesize/1000.0*pow((((va[i]+vp[i])/2.0)*vscale)/220.0,-1.0); //in Myr
        } else if (mgc[i] < 100.0) tdisso[i] = 1.0E-30;
        //printf("mass %lf r_pa %lf time %lf va %lf vp %lf\n",mgct0[i],rpa[i], tgc[i],va[i],vp[i]);
    }
    
    //Initial mass loss by SEV (BM03)
    if (gcevolu == 1) {
        for (i=0;i<nbody;i++) {
            if (mgc[i] > 1.0E-15) {
                mgc[i] *= 0.7;
            }
        }
    }
    
    
    
    /*************************
     * main integration loop *
     *************************/
    
    while (tout<tend) {
        
        
        com_correction();
        get_coeff(tout);
        get_pot(tout);
        
        if (dfmethod) get_grid();
        if (dfmethod) get_sigma();                                                           //Attention:  Velocity not exactly synchronized!
        
        if (tout>=twrite) {
            write_nbody(tout);
            twrite += DTWRITE;
        }
        
        tout+=DTOUT;
        ecorr = 0.0;
        nstepg = 0.0;
        t1 = cputim();
        
        
#pragma omp parallel for shared(mass,vx,vy,vz,x,y,z,msmbh,si_coeffs,t0,laddx,laddy,laddz,box,boy,boz,boxorb,tboxo,tboxn,mgc,rh,rh0,tscale,tdis,ncap,ngcdis,nescape,tinit,ts,vp,va,rperi,rapo,tdisso,tgc,vscale) private(a,pt,dis,accf,dt,nstep,rvo,lxo,lyo,lzo,lxn,lyn,lzn,rvn,dilaton,chit,trelax0,k,cancel_loop,xttl,vttl,tn,tv,ns0,vecr,veca,eccen) reduction(+: ecorr, nstepg)
        for (i=0;i<nbody;i++) {
            
            nstep = 0;
            
            if (mass[i] > 1.0E-30 && mgc[i] > 1.0E-30  && t0[i]<tout) {  //Modified!
                //if (mass[i]>0.0 && t0[i]<tout) {
                pt = -0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])-get_poti(i);
                cancel_loop = 0;
                
                do {
                    rvo = x[i]*vx[i]+y[i]*vy[i]+z[i]*vz[i];
                    lxo = y[i]*vz[i]-z[i]*vy[i];
                    lyo = z[i]*vx[i]-x[i]*vz[i];
                    lzo = x[i]*vy[i]-y[i]*vx[i];
                    
                    // Reduce timesteps for increased accuracy close to the SMBH
                    
                    dis = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
                    
                    //####### Integrator selector ########
                    if (rkselect == 0) {
                        // ======= Begin of MT integrator =======
                        
                        if (dis>msmbh) accf = 1.0;
                        else if (dis>msmbh/10.0) accf = 0.1;
                        else accf = 0.03;
                        
                        for (k = 0; k < si_coeffs[ORDER]; k++) {
                            dt = si_c[ORDER][k]*accf*H/(0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])+pt);
                            //if (t0[i]+dt > tout) dt = tout-t0[i];   //force all particles to be at t0 = tout at end of loop
                            t0[i] += dt;
                            x[i] += vx[i]*dt;
                            y[i] += vy[i]*dt;
                            z[i] += vz[i]*dt;
                            
                            if (si_d[ORDER][k]!=0.0) {
                                dilaton = 0.0;
                                getforce_scf(i, dilaton, a);
                                dt = -accf*H/get_poti(i)*si_d[ORDER][k];
                                vx[i] += dt*a[0];
                                vy[i] += dt*a[1];
                                vz[i] += dt*a[2];
                            }
                        }
                        
                        nstep++;
                    }
                    //======= End of MT integrator ========
                    else if (rkselect == 1) {
                        //======= TTL integrator ==============
                        xttl[0]=x[i];   xttl[1]=y[i];   xttl[2]=z[i];
                        vttl[0]=vx[i];  vttl[1]=vy[i];  vttl[2]=vz[i];
                        tv = t0[i];
                        get_ttl(i,xttl,vttl);
                        tn = t0[i];
                        nstep++;
                        x[i]=xttl[0];   y[i]=xttl[1];     z[i]=xttl[2];
                        vx[i]=vttl[0];  vy[i]=vttl[1];    vz[i]=vttl[2];
                        //======= End of TTL ==================
                    }
                    rvn = x[i]*vx[i]+y[i]*vy[i]+z[i]*vz[i];
                    lxn = y[i]*vz[i]-z[i]*vy[i];
                    lyn = z[i]*vx[i]-x[i]*vz[i];
                    lzn = x[i]*vy[i]-y[i]*vx[i];
                    laddx[i] += lxn;
                    laddy[i] += lyn;
                    laddz[i] += lzn;
                    
                    if ((lxo*lxn) < 0.0) box[i] = 1;
                    if ((lyo*lyn) < 0.0) boy[i] = 1;
                    if ((lzo*lzn) < 0.0) boz[i] = 1;
                    if ((box[i]*boy[i]*boz[i]) == 1) {
                        boxorb[i] = 1;
                        tboxo[i]   = t0[i]-tboxn[i];
                        tboxn[i]   = t0[i];
                        box[i]     = 0;
                        boy[i]     = 0;
                        boz[i]     = 0;
                    }
                    if ((boxorb[i] == 1) && ((t0[i]-tboxn[i]) > (5*tboxo[i]))) boxorb[i]=0;
                    
                    if (rvo<0.0 && rvn > 0.0 && mass[i] > 0.0 && check_smbh_capture(i)) {
                        ncap++;
                        ecorr += mass[i]*(get_poti(i)+0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]));
                        mass[i] = 0.0;
                        cancel_loop = 1;
                    }
                    
                    if (gcevolu == 1 && mgc[i] > 1.0E-30 && mass[i] > 1.0E-30) { //Mass loss through relaxation in tidal fields
                        /*chit = 3.0*pow(t0[i]*tscale/2.0,-0.3);   //Size evolution (for later publications)
                         trelax0 = 0.134998*pow(mgc[i]/0.5,0.5)*pow(rh0[i],1.5)/(0.046378875*log(0.22*mgc[i])); //in Myr
                         rh[i]  = rh0[i]*pow(pow(t0[i]*tscale/2.0,0.14)+pow(chit*t0[i]*tscale/trelax0,1.3333333),0.5); */
                        
                        if (rvo<0.0 && rvn > 0.0) {                   // Pericenter
                            rperi[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
                            vecr[0]     = x[i];
                            vecr[1]     = y[i];
                            vecr[2]     = z[i];
                            get_forcettl(vecr, veca);
                            vp[i]       = sqrt(rperi[i]*sqrt(veca[0]*veca[0]+veca[1]*veca[1]+veca[2]*veca[2]));
                            tgc[i]      = 0.0;
                            if (mgc[i] > 100.0) {
                                ns0 = mgct0[i]/0.64;
                                tdisso[i] = 1.91*pow(ns0/(log(0.02*ns0)),0.75)*rperi[i]*scalesize/1000.0*pow((((va[i]+vp[i])/2.0)*vscale)/220.0,-1.0); //in Myr
                            }
                        }
                        if (rvo>0.0 && rvn < 0.0) {                   // Apocenter
                            rapo[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
                            vecr[0]     = x[i];
                            vecr[1]     = y[i];
                            vecr[2]     = z[i];
                            get_forcettl(vecr, veca);
                            va[i]       = sqrt(rapo[i]*sqrt(veca[0]*veca[0]+veca[1]*veca[1]+veca[2]*veca[2]));
                            tgc[i]      = 0.0;
                            if (mgc[i] > 100.0) {
                                ns0 = mgct0[i]/0.64;
                                tdisso[i] = 1.91*pow(ns0/(log(0.02*ns0)),0.75)*rapo[i]*scalesize/1000.0*pow((((va[i]+vp[i])/2.0)*vscale)/220.0,-1.0); //in Myr
                            }
                        }
                        if (rapo[i] >= rperi[i] && rapo[i] > 1.0E-20) eccen = (rapo[i] - rperi[i])/(rapo[i] + rperi[i]);
                        tgc[i] += tn-tv;
                        mgc[i] -= (tn-tv)*tscale*mgct0[i]/tdisso[i];
                        if (mgc[i] <= 100.0) {                                                  //Dissolution criterium
                            ngcdis++;
                            mgc[i] = 0.0;
                            if (boxorb[i] == 1) {
                                printf("Dissolution of GC on Box orbit: %7i Disruptions: %i time: %10.5lf physical time [Myr]  %10.5lf r: %10.6lf\n",i,ngcdis,t0[i],t0[i]*tscale,sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]));
                            } else {
                                printf("Dissolution of GC on non Box orbit: %7i Disruptions: %i time: %10.5lf physical time [Myr] %10.5lf  r: %10.6lf\n",i,ngcdis,t0[i],t0[i]*tscale,sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]));
                            }
                        }
                        //printf("time %lf mass %lf Tdiss %lf eccentricity %lf\n",t0[i]*tscale, mgc[i],tdisso[i], eccen);
                        
                        
                    }
                    
                    if (tdis == 1 || tdis == 2) {
                        if (t0[i]*tscale > tinit && mgc[i] > 1.0E-30 && mass[i] > 1.0E-30 && check_gc_disruption(i)) {                // Disruption by tidal Shock
                            ngcdis++;
                            mgc[i] = 0.0;
                        }
                    }
                    
                    if (mass[i] > 0.0 && dis>RMAX) {
                        nescape++;
                        ecorr += mass[i]*(get_poti(i)+0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]));
                        mass[i] = 0.0;
                        cancel_loop = 1;
                    }
                    
                    
                    
                    
                    
                    //*******TEST******
                    /*if (mgc[i] < 0.0 && t0[i] > ts[i] && i==1) {
                     printf("time: %lf Radius: %lf GC_Mass %lf\n ",t0[i]*tscale,sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]),mgc[i]);
                     ts[i]+=1.0;
                     }*/
                    //*****TESTEND*****
                    
                    
                    
                } while (t0[i]<tout && cancel_loop == 0);
            }
            
            //     printf("thread: %2i  particle: %5i  nsteps: %10i  time: %10lf\n",omp_get_thread_num(),i,nstep,cputim());
            nstepg += nstep;
            
        }
        
        t2 = cputim();
        printf("\nTime: %lf Integration time: %lf per time unit: %lf\n",tout,t2-t1,(t2-t1)/DTOUT);
        
        nstepg = (double) nstepg/(1.0*nbody*DTOUT);
        printf("\nAverage number of steps/particle & time unit: %lf\n",nstepg);
        
        fflush(stdout);
        
        // Write restart File in case of crash
        write_restart();
        
    }
    write_nbody(tout);
    get_pot(tout);
    fclose(outputtest);
}


double get_poti(int i) {
    
    int l,m,n;
    double clm[LMAX][LMAX],dlm[LMAX][LMAX];
    double r,xi,phi,cosphi,costheta,pot;
    double phibar[LMAX][NMAX],p[LMAX][LMAX],c[LMAX][NMAX];
    double fac1,fac2;
    
    r    = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
    xi   = (r-1.0)/(r+1.0);
    cosphi = x[i]/sqrt(x[i]*x[i]+y[i]*y[i]);
    phi = acos(cosphi);
    if (y[i]<0.0) phi=2.0*PI-phi;
    costheta = z[i]/r;
    
    for (l=0;l<LMAX;l++) {
        fac1 = 2.0*l+1.5;
        c[l][0] = (double) 1.0;
        c[l][1] = (double) 2.0*fac1*xi;
        for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*(2.0*(n-1.0+fac1)*xi*c[l][n-1]-(n+2.0*fac1-2.0)*c[l][n-2]);
    }
    
    fac1 = -sqrt(4.0*PI)/(1.0+r);
    
    for (l=0;l<LMAX;l++) {
        fac2 = fac1*pow(r,l)/pow(1.0+r,2.0*l);
        for (n=0;n<NMAX;n++) phibar[l][n]=fac2*c[l][n];
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            clm[l][m] = 0.0;
            dlm[l][m] = 0.0;
        }
    }
    
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac1 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac2 = fac1*abar[l][n]*phibar[l][n];
                clm[l][m] += fac2*clmfull[l][m][n];
                dlm[l][m] += fac2*slmfull[l][m][n];
            }
        }
    }
    
    
    for (m=0;m<LMAX;m++) {
        p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(1.0-costheta*costheta,0.5*m);
        if (m<LMAX-1) p[m+1][m] = (double) costheta*(2.0*m+1.0)*p[m][m];
        for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costheta*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
    }
    
    pot = 0.0;
    
    for (l=0;l<LMAX;l++)
        for (m=0;m<=l;m++)
            pot += (double) p[l][m]*(clm[l][m]*cos(m*phi)+dlm[l][m]*sin(m*phi));
    
    //  printf("i: %i v: %lf\n",i,pot);
    //  fflush(stdout);
    
    // ADD SMBH contribution
    pot += get_phi_smbh(i);
    
    //  printf("i: %i n: %lf\n",i,pot);
    //  fflush(stdout);
    
    return(pot);
    
}


double get_densiti(int i) {
    
    int l,m,n;
    double clm[LMAX][LMAX],dlm[LMAX][LMAX];
    double r,xi,phi,cosphi,costheta,den;
    double denbar[LMAX][NMAX],p[LMAX][LMAX],c[LMAX][NMAX];
    double fac1,fac2;
    
    r    = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
    xi   = (r-1.0)/(r+1.0);
    cosphi = x[i]/sqrt(x[i]*x[i]+y[i]*y[i]);
    phi = acos(cosphi);
    if (y[i]<0.0) phi=2.0*PI-phi;
    costheta = z[i]/r;
    
    for (l=0;l<LMAX;l++) {
        fac1 = 2.0*l+1.5;
        c[l][0] = (double) 1.0;
        c[l][1] = (double) 2.0*fac1*xi;
        for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*(2.0*(n-1.0+fac1)*xi*c[l][n-1]-(n+2.0*fac1-2.0)*c[l][n-2]);
    }
    
    fac1 = 1.0/(sqrt(PI)*r*pow(1.0+r,3));
    
    for (l=0;l<LMAX;l++) {
        fac2 = fac1*pow(r,l)/pow(1.0+r,2.0*l);
        for (n=0;n<NMAX;n++) denbar[l][n]=(0.5*n*(n+4.0*l+3.0)+(l+1.0)*(2.0*l+1.0))*fac2*c[l][n];
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            clm[l][m] = 0.0;
            dlm[l][m] = 0.0;
        }
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac1 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac2 = fac1*abar[l][n]*denbar[l][n];
                clm[l][m] += fac2*clmfull[l][m][n];
                dlm[l][m] += fac2*slmfull[l][m][n];
            }
        }
    }
    
    for (m=0;m<LMAX;m++) {
        p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(1.0-costheta*costheta,0.5*m);
        if (m<LMAX-1) p[m+1][m] = (double) costheta*(2.0*m+1.0)*p[m][m];
        for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costheta*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
    }
    
    den = 0.0;
    
    for (l=0;l<LMAX;l++)
        for (m=0;m<=l;m++)
            den += (double) p[l][m]*(clm[l][m]*cos(m*phi)+dlm[l][m]*sin(m*phi));
    
    return(den);
    
}


void get_pot(double tout) {
    
    int i,l,m,n;
    double pot2;
    double clm[LMAX][LMAX],dlm[LMAX][LMAX];
    double r,xi,phi,cosphi,costheta,ekin,pot,pottot = 0.0;
    double phibar[LMAX][NMAX],p[LMAX][LMAX],c[LMAX][NMAX];
    double fac1,fac2;
    
#pragma omp parallel for shared(nbody,clmfull,slmfull,x,y,z,mass) private (fac1,fac2,clm,dlm,pot,r,xi,cosphi,phi,costheta,c,pot2,phibar,p,m,l,n) reduction (+: pottot)
    for (i=0;i<nbody;i++) {
        r    = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
        xi   = (r-1.0)/(r+1.0);
        cosphi = x[i]/sqrt(x[i]*x[i]+y[i]*y[i]);
        phi = acos(cosphi);
        if (y[i]<0.0) phi=2.0*PI-phi;
        costheta = z[i]/r;
        
        for (l=0;l<LMAX;l++) {
            fac1 = 2.0*l+1.5;
            c[l][0] = (double) 1.0;
            c[l][1] = (double) 2.0*fac1*xi;
            for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*(2.0*(n-1.0+fac1)*xi*c[l][n-1]-(n+2.0*fac1-2.0)*c[l][n-2]);
        }
        
        fac1 = -sqrt(4.0*PI)/(1.0+r);
        for (l=0;l<LMAX;l++) {
            fac2 = fac1*pow(r,l)/pow(1.0+r,2.0*l);
            for (n=0;n<NMAX;n++) phibar[l][n]=fac2*c[l][n];
        }
        
        for (l=0;l<LMAX;l++) {
            for (m=0;m<=l;m++) {
                clm[l][m] = 0.0;
                dlm[l][m] = 0.0;
            }
        }
        
        for (l=0;l<LMAX;l++) {
            for (m=0;m<=l;m++) {
                fac1 = narr[l][m];
                for (n=0;n<NMAX;n++) {
                    fac2 = fac1*abar[l][n]*phibar[l][n];
                    clm[l][m] += fac2*clmfull[l][m][n];
                    dlm[l][m] += fac2*slmfull[l][m][n];
                }
            }
        }
        
        for (m=0;m<LMAX;m++) {
            p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(1.0-costheta*costheta,0.5*m);
            if (m<LMAX-1) p[m+1][m] = (double) costheta*(2.0*m+1.0)*p[m][m];
            for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costheta*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
        }
        
        pot = 0.0;
        
        for (l=0;l<LMAX;l++)
            for (m=0;m<=l;m++)
                pot += (double) p[l][m]*(clm[l][m]*cos(m*phi)+dlm[l][m]*sin(m*phi));
        
        pot2 = 0.0;
        
        /*
         for (k=i+1;k<nbody;k++) {
         dis=sqrt(pow(x[i]-x[k],2.0)+pow(y[i]-y[k],2.0)+pow(z[i]-z[k],2.0));
         pot2 += -mass[k]/dis;
         }
         
         printf("Pot: %lf %lf \n",pot,pot2);
         
         */
        
        // ADD SMBH contribution
        phif[i]=pot+get_phi_smbh(i);
        pot += 2.0*get_phi_smbh(i);
        
        //  pot2 += get_phi_smbh(i);
        pottot += mass[i]*pot/2.0;
        //  thpot2 += mass[i]*pot2;
        
    }
    
    ekin = 0.0;
    mcom = 0.0;
    rxcom = 0.0;
    rycom = 0.0;
    rzcom = 0.0;
    vxcom = 0.0;
    vycom = 0.0;
    vzcom = 0.0;
    
    for (i=0;i<nbody;i++) {
        mcom += mass[i];
        rxcom += mass[i]*x[i];
        rycom += mass[i]*y[i];
        rzcom += mass[i]*z[i];
        vxcom += mass[i]*vx[i];
        vycom += mass[i]*vy[i];
        vzcom += mass[i]*vz[i];
        if (boxorb[i] == 1) {
            nbox++;
        }
        ekin += 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    }
    
    if (tout==0.0) etot0 = pottot+ekin;
    
    rcom=sqrt(rxcom*rxcom+rycom*rycom+rzcom*rzcom)/mcom;
    vcom=sqrt(vxcom*vxcom+vycom*vycom+vzcom*vzcom)/mcom;
    
    // printf("Energies: %lf %lf %lf\n",pottot+ekin,ecorr,etot0);
    printf("\nTime: %lf    Total energy: %12.6le    Kinetic energy: %12.6le Virial ratio 2T/|W|  %12.6le  Energy change: %10.5lf\n",tout,pottot+ekin+ecorr,ekin,-2*ekin/pottot,(pottot+ekin+ecorr-etot0)/etot0);
    printf("             Center of mass: %12.6le    CoM velocity: %12.6le   nescape: %i\n",rcom,vcom,nescape);
    if (nbox > 0) {
        printf("Total number of box orbit stars in the galaxy: %i, fraction: %lf\n",nbox,(double)nbox/(double)nbody);
        nbox=0;
    }
    
    fflush(stdout);
    
}


void get_coeff(double tout) {
    
    int i,l,m,n;
    double deltat,r,xi,phik,cosphik,costhetak;
    double phibar[LMAX][NMAX],p[LMAX][LMAX],c[LMAX][NMAX];
    double fac1,fac2;
    double *x0,*y0,*z0;
    
    x0=malloc(BMAX*sizeof(double));
    y0=malloc(BMAX*sizeof(double));
    z0=malloc(BMAX*sizeof(double));
    
    for (i=0;i<nbody;i++) {
        deltat=tout-t0[i];
        x0[i]=x[i]+deltat*vx[i];
        y0[i]=y[i]+deltat*vy[i];
        z0[i]=z[i]+deltat*vz[i];
    }
    
    //ist hier Parallelisierung nicht effektiver?
#pragma omp parallel for shared(clmfull, slmfull) private(m, n)
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++)  {
            //Parallelisierung hier macht doch wenig Sinn!?
            //#pragma omp parallel for shared (l,m,clmfull, slmfull)
            for (n=0;n<NMAX;n++) {
                clmfull[l][m][n] = 0.0;
                slmfull[l][m][n] = 0.0;
            }
        }
    }
    
    // #pragma omp parallel for shared(nbody,x,y,z,mass,clmfull,slmfull) private (fac1,fac2,l,m,n,c,phibar,r,xi,cosphik,phik,costhetak,p)
    
    for (i=0;i<nbody;i++) {
        r    = sqrt(x0[i]*x0[i]+y0[i]*y0[i]+z0[i]*z0[i]);
        xi   = (r-1.0)/(r+1.0);
        cosphik = x0[i]/sqrt(x0[i]*x0[i]+y0[i]*y0[i]);
        phik = acos(cosphik);
        if (y0[i]<0.0) phik=2.0*PI-phik;
        costhetak = z0[i]/r;
        
        for (l=0;l<LMAX;l++) {
            c[l][0] = (double) 1.0;
            c[l][1] = (double) 2.0*(2.0*l+1.5)*xi;
            fac1 = 2.0*l+1.5;
            for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*(2.0*(n-1.0+fac1)*xi*c[l][n-1]-(n+2.0*fac1-2.0)*c[l][n-2]);
        }
        
        fac1 = -sqrt(4.0*PI)/(1.0+r);
        
        for (l=0;l<LMAX;l++) {
            fac2 = fac1*pow(r,l)/pow(1.0+r,2.0*l);
            for (n=0;n<NMAX;n++) phibar[l][n]=fac2*c[l][n];
        }
        
        fac1 = sqrt(1.0-costhetak*costhetak);
        
        for (m=0;m<LMAX;m++) {
            p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(fac1,m);
            if (m<LMAX-1) p[m+1][m] = (double)  costhetak*(2.0*m+1.0)*p[m][m];
            for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costhetak*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
        }
        
#pragma omp parallel for shared (i,phik,clmfull, slmfull,mass,p,phibar) private (m,n,fac1,fac2)
        for (l=0;l<LMAX;l++) {
            for (m=0;m<=l;m++)  {
                fac1 = mass[i]*p[l][m]*cos(m*phik);
                fac2 = mass[i]*p[l][m]*sin(m*phik);
                for (n=0;n<NMAX;n++) {
                    clmfull[l][m][n] += phibar[l][n]*fac1;
                    slmfull[l][m][n] += phibar[l][n]*fac2;
                }
            }
        }
    }
    
    free(x0);
    free(y0);
    free(z0);
    
}


void getforce_direct(int par, double xi[3], double *a) {
    
    int i,k;
    double dis,dis3,epsilon=1.e-4;
    
    for (k=0;k<3;k++) a[k] = 0.0;
    
    for (i=0;i<nbody;i++) {
        dis=sqrt(pow(xi[0]-x[i],2.0)+pow(xi[1]-y[i],2.0)+pow(xi[2]-z[i],2.0)+epsilon);
        dis3=dis*dis*dis;
        if (i != par) {
            a[0] += -mass[i]*(xi[0]-x[i])/dis3;
            a[1] += -mass[i]*(xi[1]-y[i])/dis3;
            a[2] += -mass[i]*(xi[2]-z[i])/dis3;
        }
    }
    
    // ADD SMBH contribution
    get_force_smbh(xi, a);
    
}


void getforce_scf(int i, double dilaton, double *a) {
    
    int n,l,m;
    double phibar[LMAX][NMAX],dphibar[LMAX][NMAX],p[LMAX][LMAX],dp[LMAX][LMAX],c[LMAX][NMAX],cmod[LMAX][NMAX];
    double elm[LMAX][LMAX],flm[LMAX][LMAX],clm[LMAX][LMAX],dlm[LMAX][LMAX];
    double r,xi,cosphi,phi,costheta,theta;
    double ar,atheta,aphi;
    double fac1,fac2,fac3,fac4;
    double xp[3],xgc[3],xdf[3],vdf[3];
    
    xgc[0] = x[i];
    xgc[1] = y[i];
    xgc[2] = z[i];
    
    if (dilaton > 0.0 && mgc[i] > 0.0 && tdis == 2) {
        xgc[0] = dilaton*x[i];
        xgc[1] = dilaton*y[i];
        xgc[2] = dilaton*z[i];
    }
    
    r    = sqrt(xgc[0]*xgc[0]+xgc[1]*xgc[1]+xgc[2]*xgc[2]);
    xi   = (r-1.0)/(r+1.0);
    cosphi = xgc[0]/sqrt(xgc[0]*xgc[0]+xgc[1]*xgc[1]);
    phi = acos(cosphi);
    if (xgc[1]<0.0) phi=2.0*PI-phi;
    costheta = xgc[2]/r;
    theta = acos(costheta);
    
    for (l=0;l<LMAX;l++) {
        fac1 = 2.0*(2.0*l+1.5);
        c[l][0] = (double) 1.0;
        c[l][1] = (double) fac1*xi;
        for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*((2.0*n-2.0+fac1)*xi*c[l][n-1]-(n+fac1-2.0)*c[l][n-2]);
    }
    
    for (l=0;l<LMAX;l++) {
        fac1 = (double) -sqrt(4.0*PI)*pow(r,l)/pow(1.0+r,2.0*l+1.0);
        for (n=0;n<NMAX;n++) phibar[l][n]=fac1*c[l][n]*abar[l][n];
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            elm[l][m] = 0.0;
            flm[l][m] = 0.0;
            clm[l][m] = 0.0;
            dlm[l][m] = 0.0;
        }
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac3 = 0.0;
            fac4 = 0.0;
            fac2 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac1 = phibar[l][n];
                fac3 += fac1*clmfull[l][m][n];
                fac4 += fac1*slmfull[l][m][n];
            }
            clm[l][m] += fac3*fac2;
            dlm[l][m] += fac4*fac2;
        }
    }
    
    for (l=0;l<LMAX;l++) {
        cmod[l][0] = 0.0;
        fac1 = 1.0/((4.0*l+3.0)*(1.0-xi*xi));
        fac2 = 4.0*l+2.0;
        for (n=1;n<NMAX;n++) cmod[l][n]= (double) fac1*((fac2+n)*c[l][n-1]/c[l][n]-n*xi);   // from Wolfram Math World
    }
    
    
    
    for (l=0;l<LMAX;l++) {
        fac1 = (double) l/r-(2.0*l+1.0)/(1.0+r);
        fac2 = 4.0*(2.0*l+1.5)/pow(1.0+r,2.0);
        for (n=0;n<NMAX;n++) dphibar[l][n] = phibar[l][n]*(fac1+fac2*cmod[l][n]);
    }
    
    
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac3 = 0.0;
            fac4 = 0.0;
            fac2 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac1 = dphibar[l][n];
                fac3 += fac1*clmfull[l][m][n];
                fac4 += fac1*slmfull[l][m][n];
            }
            elm[l][m] += fac3*fac2;
            flm[l][m] += fac4*fac2;
        }
    }
    
    fac1 = sqrt(1.0-costheta*costheta);
    
    for (m=0;m<LMAX;m++) {
        p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(fac1,m);
        if (m<LMAX-1) p[m+1][m] = (double) costheta*(2.0*m+1.0)*p[m][m];
        for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costheta*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
    }
    
    fac1 = costheta/sin(theta);
    
    for (l=0;l<LMAX;l++) {
        dp[l][l] = (double) l*fac1*p[l][l];
        for (m=0;m<l;m++) dp[l][m] = (double) m*fac1*p[l][m]+p[l][m+1];
    }
    
    ar = 0.0; aphi = 0.0; atheta = 0.0;
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac1 = cos(m*phi);
            fac2 = sin(m*phi);
            ar     += -p[l][m]*(elm[l][m]*fac1+flm[l][m]*fac2);
            atheta += -dp[l][m]/r*(clm[l][m]*fac1+dlm[l][m]*fac2);
            aphi   += -m*p[l][m]/r/sin(theta)*(dlm[l][m]*fac1-clm[l][m]*fac2);
        }
    }
    
    a[0] = sin(theta)*cos(phi)*ar+cos(theta)*cos(phi)*atheta-sin(phi)*aphi;
    a[1] = sin(theta)*sin(phi)*ar+cos(theta)*sin(phi)*atheta+cos(phi)*aphi;
    a[2] = cos(theta)*ar-sin(theta)*atheta;
    
    //****Test for pure Hernquist Model****
    /*a[0] = -1.0*pow(0.414213562+r,-2.0)*xgc[0]/r;
     a[1] = -1.0*pow(0.414213562+r,-2.0)*xgc[1]/r;
     a[2] = -1.0*pow(0.414213562+r,-2.0)*xgc[2]/r;*/
    //*************************************
    
    //End Test
    // ADD SMBH contribution
    xp[0] = xgc[0];  xp[1] = xgc[1]; xp[2] = xgc[2];
    get_force_smbh(xp, a);
    
    // Add DYNAMICAL FRICTION // Attention: mgc[i] has physical dimension
    if (dfmethod > 0 && mgc[i] > 0.0 && dilaton == 0.0 && t0[i]*tscale > tinit) {
        xdf[0] = xgc[0], xdf[1] = xgc[1]; xdf[2] = xgc[2];
        vdf[0] = vx[i], vdf[1] = vy[i],  vdf[2] = vz[i];
        get_force_df(i, xdf, vdf, a);
    }
}


double factrl(int n) {
    
    static int ntop=4;
    static double a[33]={1.0,1.0,2.0,6.0,24.0};
    int j;
    
    if (n < 0) {printf("Negative factorial in routine factrl: %i\n",n); exit(-1);}
    if (n > 32) return exp(gammln(n+1.0));
    
    while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
    }
    
    return a[n];
    
}


double gammln(double xx) {
    
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    
    return -tmp+log(2.5066282746310005*ser/x);
    
}


void wtime() {
    
    struct timeval time1;
    int ierr=0;
    
    ierr = gettimeofday(&time1, NULL) ;
    if (ierr != 0 ) printf("bad return of gettimeofday, ierr = %d \n",ierr);
    tim[0] = time1.tv_sec;
    tim[1] = time1.tv_usec;
    
}


double cputim() {
    /*
     This timer returns the difference between current
     wall clock time and the one returned by wtime_ and
     returns a double precision number
     */
    
    struct timeval time1;
    int ierr=0;
    
    if (tim[0]==0.0 && tim[1]==0.0) wtime();
    ierr = gettimeofday(&time1, NULL);
    if (ierr != 0 ) printf("bad return of gettime of day, ierr = %d \n", ierr);
    
    return (double) (time1.tv_sec - tim[0])/60.0 + (time1.tv_usec - tim[1])/6.e7;
    
}


double get_phi_smbh(int i) {
    
    double r,phi_smbh;
    
    r = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
    phi_smbh = -msmbh/r;
    
    return phi_smbh;
    
}


void get_force_smbh(double *xi, double *a) {
    
    double dis;
    
    dis = sqrt(xi[0]*xi[0]+xi[1]*xi[1]+xi[2]*xi[2]);
    a[0] -= msmbh/(dis*dis*dis)*xi[0];
    a[1] -= msmbh/(dis*dis*dis)*xi[1];
    a[2] -= msmbh/(dis*dis*dis)*xi[2];
    
}


void get_grid() {
    
    double l;
    int i, j;
    for (i=0;i<nbody;i++) grid_cell_nr[i] = 0;
    
    
    for (j=1; j<=GRIDDEPTH; j++) {  //Attention, j runs from 1 to GRIDDEPTH
        
        l = pow(0.6,j-1)*RMAX/(1.0*GRIDCELLS); //half grid-cell length is half simulation-volume length divided by GRIDCELLS
        
        //printf("\nHalf box length = %f\n",l);
        
#pragma omp parallel for shared(j, l, grid_cell_nr, x, y, z)
        for(i=0;i<nbody;i++) {
            if ((!grid_cell_nr[i]) && (mass[i]>0.0)) {
                if ((x[i]<l) && (x[i]>-l)) {
                    if ((y[i]<l) && (y[i]>-l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 0; //0
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-1; //98
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-2; //97
                        }
                    } else if ((y[i]>=l) && (y[i]<3.0*l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-3; //96
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-4; //95
                        }
                    } else if ((y[i]>=3.0*l) && (y[i]<=5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-5; //94
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-6; //93
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-7; //92
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-8; //91
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-9; //90
                        }
                    } else if ((y[i]<=-l) && (y[i]>-3.0*l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-10; //89
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-11; //88
                        }
                    } else if ((y[i]<=-3.0*l) && (y[i]>=-5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-12; //87
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-13; //86
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-14; //85
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-15; //84
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-16; //83
                        }
                    }
                } else if ((x[i] >= l) && (x[i] < 3.0*l)) {
                    if ((y[i]<l) && (y[i]>-l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-17; //82
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-18; //81
                        }
                    } else if ((y[i]>=l) && (y[i]<3.0*l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-19; //80
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-20; //79
                        }
                    } else if ((y[i]>=3.0*l) && (y[i]<=5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-21; //78
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-22; //77
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-23; //76
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-24; //75
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-25; //74
                        }
                    } else if ((y[i]<=-l) && (y[i]>-3.0*l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-26; //73
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-27; //72
                        }
                    } else if ((y[i]<=-3.0*l) && (y[i]>=-5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-28; //71
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-29; //70
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-30; //69
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-31; //68
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-32; //67
                        }
                    }
                } else if ((x[i] >= 3.0*l) && (x[i] <= 5.0*l)) {
                    if ((y[i]<l) && (y[i]>-l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-33; //66
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-34; //65
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-35; //64
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-36; //63
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-37; //62
                        }
                    } else if ((y[i]>=l) && (y[i]<3.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-38; //61
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-39; //60
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-40; //59
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-41; //58
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-42; //57
                        }
                    } else if ((y[i]>=3.0*l) && (y[i]<=5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-43; //56
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-44; //55
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-45; //54
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-46; //53
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-47; //52
                        }
                    } else if ((y[i]<=-l) && (y[i]>-3.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-48; //51
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-49; //50
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-50; //49
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-51; //48
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-52; //47
                        }
                    } else if ((y[i]<=-3.0*l) && (y[i]>=-5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-53; //46
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-54; //45
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-55; //44
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-56; //43
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-57; //42
                        }
                    }
                } else if ((x[i] <= -l) && (x[i] > -3.0*l)) {
                    if ((y[i]<l) && (y[i]>-l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-58; //41
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-59; //40
                        }
                    } else if ((y[i]>=l) && (y[i]<3.0*l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-60; //39
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-61; //38
                        }
                    } else if ((y[i]>=3.0*l) && (y[i]<=5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-62; //37
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-63; //36
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-64; //35
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-65; //34
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-66; //33
                        }
                    } else if ((y[i]<=-l) && (y[i]>-3.0*l)) {
                        if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-67; //32
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-68; //31
                        }
                    } else if ((y[i]<=-3.0*l) && (y[i]>=-5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-69; //30
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-70; //29
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-71; //28
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-72; //27
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-73; //26
                        }
                    }
                } else if ((x[i] <= -3.0*l) && (x[i] >= -5.0*l)) {
                    if ((y[i]<l) && (y[i]>-l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-74; //25
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-75; //24
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-76; //23
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-77; //22
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-78; //21
                        }
                    } else if ((y[i]>=l) && (y[i]<3.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-79; //20
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-80; //19
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-81; //18
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-82; //17
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-83; //16
                        }
                    } else if ((y[i]>=3.0*l) && (y[i]<=5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-84; //15
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-85; //14
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-86; //13
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-87; //12
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-88; //11
                        }
                    } else if ((y[i]<=-l) && (y[i]>-3.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-89; //10
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-90; //9
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-91; //8
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-92; //7
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-93; //6
                        }
                    } else if ((y[i]<=-3.0*l) && (y[i]>=-5.0*l)) {
                        if ((z[i]<l) && (z[i]>-l)) {
                            grid_cell_nr[i] = 98*j+1-94; //5
                        } else if ((z[i]>=l) && (z[i]<3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-95; //4
                        } else if ((z[i]>=3.0*l) && (z[i]<=5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-96; //3
                        } else if ((z[i]<=-l) && (z[i]>-3.0*l)) {
                            grid_cell_nr[i] = 98*j+1-97; //2
                        } else if ((z[i]<=-3.0*l) && (z[i]>=-5.0*l)) {
                            grid_cell_nr[i] = 98*j+1-98; //1
                        }
                    }
                }
            }
        }
    }
    
    for (i=0;i<98*GRIDDEPTH+1;i++) grid_cell_occupation[i] = 0;
    for (i=0;i<nbody;i++) grid_cell_occupation[grid_cell_nr[i]] ++;
    
    for (i=0;i<98*GRIDDEPTH+1;i++) {
        grid_coordinates[i][0] = 0.0;//x
        grid_coordinates[i][1] = 0.0;//y
        grid_coordinates[i][2] = 0.0;//z
        grid_coordinates[i][3] = 0.0;//M
    }
    for (i=0;i<nbody;i++) {
        if (mass[i]>0.0) {
            grid_coordinates[grid_cell_nr[i]][0] += mass[i]*x[i];
            grid_coordinates[grid_cell_nr[i]][1] += mass[i]*y[i];
            grid_coordinates[grid_cell_nr[i]][2] += mass[i]*z[i];
            grid_coordinates[grid_cell_nr[i]][3] += mass[i];
        }
    }
    for (i=0;i<98*GRIDDEPTH+1;i++) {
        if (grid_cell_occupation[i] > 2) {
            grid_coordinates[i][0] /= 1.0*grid_coordinates[i][3];//mx
            grid_coordinates[i][1] /= 1.0*grid_coordinates[i][3];//my
            grid_coordinates[i][2] /= 1.0*grid_coordinates[i][3];//mz
        } else {
            grid_coordinates[i][0] = 0.0;
            grid_coordinates[i][1] = 0.0;
            grid_coordinates[i][2] = 0.0;
            grid_coordinates[i][3] = 0.0;
        }
    }
    
    
    /*	int tot_occupation=0;
     for (i=0;i<98*GRIDDEPTH+1;i++) {
     printf ("%i\t%i\t\t%f\t%f\t%f\t%f\n", i, grid_cell_occupation[i], grid_coordinates[i][0], grid_coordinates[i][1], grid_coordinates[i][2], grid_coordinates[i][3]);
     tot_occupation += grid_cell_occupation[i];
     }
     printf("\nTotal occupation: %i\n", tot_occupation);
     */
    return;
    
}


void jacobi(float **am, int n, float d[], float **v, int *nrot) {
    
    int j,iq,ip,i;
    float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
    
    b=vector(1,n);
    z=vector(1,n);
    
    for (ip=1;ip<=n;ip++) {
        for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    }
    
    for (ip=1;ip<=n;ip++) {
        b[ip]=d[ip]=am[ip][ip];
        z[ip]=0.0;
    }
    
    *nrot=0;
    
    for (i=1;i<=50;i++) {
        sm=0.0;
        for (ip=1;ip<=n-1;ip++) {
            for (iq=ip+1;iq<=n;iq++)
                sm += fabs(am[ip][iq]);
        }
        if (sm == 0.0) {
            free_vector(z,1,n);
            free_vector(b,1,n);
            return;
        }
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
            tresh=0.0;
        for (ip=1;ip<=n-1;ip++) {
            for (iq=ip+1;iq<=n;iq++) {
                g=100.0*fabs(am[ip][iq]);
                if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) am[ip][iq]=0.0;
                else if (fabs(am[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if ((float)(fabs(h)+g) == (float)fabs(h)) t=(am[ip][iq])/h;
                    else {
                        theta=0.5*h/(am[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*am[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    am[ip][iq]=0.0;
                    
                    for (j=1;j<=ip-1;j++) {
                        ROTATE(am,j,ip,j,iq)
                    }
                    
                    for (j=ip+1;j<=iq-1;j++) {
                        ROTATE(am,ip,j,j,iq)
                    }
                    
                    for (j=iq+1;j<=n;j++) {
                        ROTATE(am,ip,j,iq,j)
                    }
                    
                    for (j=1;j<=n;j++) {
                        ROTATE(v,j,ip,j,iq)
                    }
                    
                    ++(*nrot);
                    
                }
            }
        }
        
        for (ip=1;ip<=n;ip++) {
            b[ip] += z[ip];
            d[ip]=b[ip];
            z[ip]=0.0;
        }
    }
    
    nrerror("Too many iterations in routine jacobi");
    
}


void eigsrt(float d[], float **v, int n) {
    
    int k,j,i;
    float p;
    for (i=1;i<n;i++) {
        p=d[k=i];
        for (j=i+1;j<=n;j++) if (d[j] >= p) p=d[k=j];
        if (k != i) {
            d[k]=d[i];
            d[i]=p;
            for (j=1;j<=n;j++) {
                p=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=p;
            }
        }
    }
    
}


void get_sigma() {
    
    int i, k, celltot, nrot;
    celltot = (98*GRIDDEPTH+1);
    double v1m[celltot], v2m[celltot], v3m[celltot];
    double v1v1m[celltot],v1v2m[celltot],v1v3m[celltot],v2v2m[celltot],v2v3m[celltot],v3v3m[celltot];
    double count[celltot];
    
    float **am_r;
    am_r = (float **)calloc(4, sizeof(float *));
    for (i=0;i<4;i++) am_r[i] = (float *)calloc(4,sizeof(float));
    
    float **v;
    v = (float **)calloc(4, sizeof(float *));
    for (i=0;i<4;i++) v[i] = (float *)calloc(4,sizeof(float));
    
    float d[4];
    
    
    for (k=0;k<celltot;k++) {   //initialize grid cell values
        
        count[k] = 0.0; //avoid division by zero
        v1m[k]=0.0;
        v2m[k]=0.0;
        v3m[k]=0.0;
        v1v1m[k]=0.0;
        v1v2m[k]=0.0;
        v1v3m[k]=0.0;
        v2v2m[k]=0.0;
        v2v3m[k]=0.0;
        v3v3m[k]=0.0;
        sigma_1[k]=0.0;
        sigma_2[k]=0.0;
        sigma_3[k]=0.0;
        
    }
    
    
    for (i=0;i<nbody;i++) {
        if (mass[i] > 0.0) {
            count[grid_cell_nr[i]] += 1.0;
            v1m[grid_cell_nr[i]] += vx[i];
            v2m[grid_cell_nr[i]] += vy[i];
            v3m[grid_cell_nr[i]] += vz[i];
            v1v1m[grid_cell_nr[i]] += vx[i]*vx[i];
            v1v2m[grid_cell_nr[i]] += vx[i]*vy[i];
            v1v3m[grid_cell_nr[i]] += vx[i]*vz[i];
            v2v2m[grid_cell_nr[i]] += vy[i]*vy[i];
            v2v3m[grid_cell_nr[i]] += vy[i]*vz[i];
            v3v3m[grid_cell_nr[i]] += vz[i]*vz[i];
        }
    }
    
    
    for (k=0;k<celltot;k++) {
        
        if (count[k]>2) {
            am_r[1][1] = (v1v1m[k]-v1m[k]/count[k]*v1m[k])/count[k];
            am_r[1][2] = (v1v2m[k]-v1m[k]/count[k]*v2m[k])/count[k];
            am_r[1][3] = (v1v3m[k]-v1m[k]/count[k]*v3m[k])/count[k];
            am_r[2][1] = am_r[1][2];
            am_r[2][2] = (v2v2m[k]-v2m[k]/count[k]*v2m[k])/count[k];
            am_r[2][3] = (v2v3m[k]-v2m[k]/count[k]*v3m[k])/count[k];
            am_r[3][1] = am_r[1][3];
            am_r[3][2] = am_r[2][3];
            am_r[3][3] = (v3v3m[k]-v3m[k]/count[k]*v3m[k])/count[k];
            
            jacobi(am_r, 3, d, v, &nrot);
            eigsrt(d, v, 3);
            
            if (d[1]>=0.0) sigma_1[k] = sqrt(d[1]);
            else sigma_1[k] = 0.0;
            
            if (d[2]>=0.0) sigma_2[k] = sqrt(d[2]);
            else sigma_2[k] = 0.0;
            
            if (d[3]>=0.0) sigma_3[k] = sqrt(d[3]);
            else sigma_3[k] = 0.0;
            
            alpha_euler[k] =  atan2(v[1][1],-v[2][1]);   // Pruefe die Eulerwinkel besser nochmal nach
            beta_euler[k] = acos(v[3][1]);
            gamma_euler[k] = atan2(v[3][2], v[3][3]);
            for (i=0;i<3;i++) {
                e1[k][i] = v[i+1][1];
                e2[k][i] = v[i+1][2];
                e3[k][i] = v[i+1][3];
            }
            
        } else {
            sigma_1[k] = 0.0;
            sigma_2[k] = 0.0;
            sigma_3[k] = 0.0;
            e1[k][0] = 1.0;
            e2[k][0] = 0.0;
            e3[k][0] = 0.0;
            e1[k][1] = 0.0;
            e2[k][1] = 1.0;
            e3[k][1] = 0.0;
            e1[k][2] = 0.0;
            e2[k][2] = 0.0;
            e3[k][2] = 1.0;
        }
        //		printf("Eigenwerte[%i] :", k);
        //		printf("%lf %lf %lf\n\n", d[1],d[2],d[3]);
        //		printf("Eigenvektoren[%i]:",k);
        //		printf("%lf %lf %lf\n", e1[k][0],e1[k][1],e1[k][2]);
        //		printf("%lf %lf %lf\n", e2[k][0],e2[k][1],e2[k][2]);
        //		printf("%lf %lf %lf\n\n", e3[k][0],e3[k][1],e3[k][2]);
        
    }
    
    for (i=0;i<4;i++) free(am_r[i]);
    free(am_r);
    
    for (i=0;i<4;i++) free(v[i]);
    free(v);
    
    //	printf ("Cell Nr\tsigma_1\tsigma_2\tsigma_3\t\talpha\tbeta\tgamma\n");
    //	for (i=0;i<98*GRIDDEPTH+1;i++) printf ("%i\t%f\t%f\t%f\t\t%f\t%f\t%f\n", i, sigma_1[i], sigma_2[i],sigma_3[i],alpha_euler[i],beta_euler[i],gamma_euler[i]);
    
    
}


void get_force_df(int i, double *xdf, double *vdf, double *a) {
    
    int k;
    double dis,vtilde[3], costheta[3];
    double adf[3] = {0.0,0.0,0.0};
    double gamma[3], gammac, vgc,lamb,sigmai,chiparam,bracketp;
    
    dis = sqrt(xdf[0]*xdf[0]+xdf[1]*xdf[1]+xdf[2]*xdf[2]);
    
    vgc = sqrt(vdf[0]*vdf[0]+vdf[1]*vdf[1]+vdf[2]*vdf[2]);
    
    
    if (rh[i] > (mgc[i]*scalesize/scalemass)) {
        lamb = log(dis*scalesize/rh[i]);
    } else {
        lamb = log(dis*scalemass/mgc[i]);
    }
    if (lamb < 0.0) lamb = 0.0;
    
    
    
    double rtemp;
    double sigma_1_inter = 0.0;
    double sigma_2_inter = 0.0;
    double sigma_3_inter = 0.0;
    double e1_inter[3] = {0.0, 0.0, 0.0};
    double e2_inter[3] = {0.0, 0.0, 0.0};
    double e3_inter[3] = {0.0, 0.0, 0.0};
    double weight = 0.0;
    
    for (k=0;k<98*GRIDDEPTH+1;k++) {
        if (grid_cell_occupation[k]>2) {
            //squared radial distance of GC from grid-box-COM
            rtemp = sqrt(pow(xdf[0]-grid_coordinates[k][0],2)+pow(xdf[1]-grid_coordinates[k][1],2)+pow(xdf[2]-grid_coordinates[k][2],2));
            if (rtemp < dis) {
                rtemp = pow(rtemp,64.0); //64 orginal
                sigma_1_inter += sigma_1[k]/rtemp;
                sigma_2_inter += sigma_2[k]/rtemp;
                sigma_3_inter += sigma_3[k]/rtemp;
                e1_inter[0] += e1[k][0]/rtemp;
                e1_inter[1] += e1[k][1]/rtemp;
                e1_inter[2] += e1[k][2]/rtemp;
                e2_inter[0] += e2[k][0]/rtemp;
                e2_inter[1] += e2[k][1]/rtemp;
                e2_inter[2] += e2[k][2]/rtemp;
                e3_inter[0] += e3[k][0]/rtemp;
                e3_inter[1] += e3[k][1]/rtemp;
                e3_inter[2] += e3[k][2]/rtemp;
                weight += 1.0/rtemp;
            }
        }
    }
    if (weight>0.0) {
        sigma_1_inter /= 1.0*weight;
        sigma_2_inter /= 1.0*weight;
        sigma_3_inter /= 1.0*weight;
        e1_inter[0] /= 1.0*weight;
        e1_inter[1] /= 1.0*weight;
        e1_inter[2] /= 1.0*weight;
        e2_inter[0] /= 1.0*weight;
        e2_inter[1] /= 1.0*weight;
        e2_inter[2] /= 1.0*weight;
        e3_inter[0] /= 1.0*weight;
        e3_inter[1] /= 1.0*weight;
        e3_inter[2] /= 1.0*weight;
    } else {
        sigma_1_inter = 1.0E10;
        sigma_2_inter = 1.0E10;
        sigma_3_inter = 1.0E10;
        e1_inter[0] = 1.0;
        e1_inter[1] = 0.0;
        e1_inter[2] = 0.0;
        e2_inter[0] = 0.0;
        e2_inter[1] = 1.0;
        e2_inter[2] = 0.0;
        e3_inter[0] = 0.0;
        e3_inter[1] = 0.0;
        e3_inter[2] = 1.0;
    }
    
    
    if (dfmethod == 1) {
        // standard Chandrasekhar
        
        sigmai = (1.0/sqrt(3.0))*sqrt(sigma_1_inter*sigma_1_inter+
                                      sigma_2_inter*sigma_2_inter+
                                      sigma_3_inter*sigma_3_inter);
        chiparam = vgc/(1.414213562*sigmai);
        gammac = 12.566370616*lamb*mgc[i]*get_densiti(i)/(scalemass*pow(vgc,3));
        bracketp = erf(chiparam)-2.0*chiparam/1.772453851*exp(-(chiparam*chiparam));
        for (k=0;k<3;k++) {
            adf[k] = -gammac*bracketp*vdf[k];
            //	  adf[k] = -0.005524446*vdf[k]/(pow(dis,2));
        }
        
    } else if (dfmethod == 2) {
        // triaxial generalization of Chandrasekhar
        
        //transform velocity vector to eigensystem of corresponding cell
        
        costheta[0] = (e1_inter[0]*vdf[0]+e1_inter[1]*vdf[1]+e1_inter[2]*vdf[2])/(sqrt(e1_inter[0]*e1_inter[0]+e1_inter[1]*e1_inter[1]+e1_inter[2]*e1_inter[2])*vgc);
        costheta[1] = (e2_inter[0]*vdf[0]+e2_inter[1]*vdf[1]+e2_inter[2]*vdf[2])/(sqrt(e2_inter[0]*e2_inter[0]+e2_inter[1]*e2_inter[1]+e2_inter[2]*e2_inter[2])*vgc);
        costheta[2] = (e3_inter[0]*vdf[0]+e3_inter[1]*vdf[1]+e3_inter[2]*vdf[2])/(sqrt(e3_inter[0]*e3_inter[0]+e3_inter[1]*e3_inter[1]+e3_inter[2]*e3_inter[2])*vgc);
        
        vtilde[0] = vgc*costheta[0];
        vtilde[1] = vgc*costheta[1];
        vtilde[2] = vgc*costheta[2];
        
        double sigma2[3];
        
        sigma2[0] = sigma_1_inter*sigma_1_inter;
        sigma2[1] = sigma_2_inter*sigma_2_inter;
        sigma2[2] = sigma_3_inter*sigma_3_inter;
        for (k=0;k<3;k++) {
            
            gammac = 5.0132565493*lamb*mgc[i]*get_densiti(i)/(scalemass*pow(sigma_1_inter,3));
            gamma[k] = gammac*get_integral(i, k, vtilde, sigma2); //Attention: sigma_1>sigma_2>sigma_3 as well as eigenvalues e1>e2>e3
        }
        //printf("\ngammas: %f\t%f\t%f\n", gamma[0], gamma[1], gamma[2]);
        
        for (k=0;k<3;k++) {
            //			adf[k] = -gamma[0]*vtilde[0]*e1[grid_cell_nr[i]][k]-gamma[1]*vtilde[1]*e2[grid_cell_nr[i]][k]-gamma[2]*vtilde[2]*e3[grid_cell_nr[i]][k];
            adf[k] = -gamma[0]*vtilde[0]*e1_inter[k]-gamma[1]*vtilde[1]*e2_inter[k]-gamma[2]*vtilde[2]*e3_inter[k];
        }
        
    }
    
    if (sqrt(adf[0]*adf[0]+adf[1]*adf[1]+adf[2]*adf[2])>sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])) {
        printf("WARNING: |F_DF| > |F_GRAV|! %f\t%f\t%f\t\t%f\t%f\t%f\n", adf[0], adf[1], adf[2], a[0], a[1], a[2]);
        printf("RESCALING F_DF to F_GRAV\n");
        for (i=0;i<3;i++) adf[i] /= sqrt(adf[0]*adf[0]+adf[1]*adf[1]+adf[2]*adf[2])/sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); //saftey measure for making sure that DF force is always <= gravitational force
    }
    
    a[0] += adf[0];
    a[1] += adf[1];
    a[2] += adf[2];
    
}


double get_integral(int i, int k, double *vtilde, double *sigma2) {
    
    int j,l;
    double s;
   	double eps2[3], vtilde2[3]; ///Original
    for (j=0;j<3;j++) eps2[j] = sigma2[j]/sigma2[0];
    for (j=0;j<3;j++) vtilde2[j] = vtilde[j]*vtilde[j];
    
    
    
    double p[10];
    for (j=0;j<3;j++) p[j] = sigma2[j];
    for (j=0;j<3;j++) p[3+j] = eps2[j];//eps[j]^2
    for (j=0;j<3;j++) p[6+j] = vtilde[j]*vtilde[j];//vtilde[j]^2
    
    p[9] = sigma2[k]/sigma2[0];
    
    
    
    s=0.0;
    double ug=0.0;
    double og=100.0;
    int NSTEPSI = 200;
    double br;
    double gr;
    br=(og-ug)/NSTEPSI;
    for (l=1;l<=NSTEPSI;l++) {
        gr=ug+(l-1)*br;
        s+=qgaus(logkernel,gr,gr+br,p);
        
    }
    //  printf("integral %f\n",s);
    //	printf("value integral = %f\n",s);
    
    return s;
    
}


double logkernel(double lu, double *p) {
    
    double u,s;
    s=exp(lu);
    u=s-1.0;
    return (exp(-0.5*p[6]/p[0]/(1.0+u)-0.5*p[7]/p[1]/(p[4]+u)-0.5*p[8]/p[2]/(p[5]+u))/((p[9]+u)*sqrt((1.0+u)*(p[4]+u)*(p[5]+u))))*s;
}



double qgaus(double (*func)(double, double*), double a, double b, double *p)
{
    int j;
    double xr,xm,dx,s;
    static double x[]={0.0,0.1488743389,0.4333953941,
        0.6794095682,0.8650633666,0.9739065285};
    static double w[]={0.0,0.2955242247,0.2692667193,
        0.2190863625,0.1494513491,0.0666713443};
    
    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;
    for (j=1;j<=5;j++) {
        dx=xr*x[j];
        s += w[j]*((*func)(xm+dx,p)+(*func)(xm-dx,p));
    }
    return s *= xr;
}





int check_smbh_capture(int i) {
    
    double r,etot,ltot,lmx,lmy,lmz,x1,x2,rmin;
    
    r = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
    etot=-msmbh/r+0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    lmx = (y[i]*vz[i]-z[i]*vy[i]);
    lmy = (z[i]*vx[i]-x[i]*vz[i]);
    lmz = (x[i]*vy[i]-y[i]*vx[i]);
    ltot = sqrt(lmx*lmx+lmy*lmy+lmz*lmz);
    x1 = -msmbh/(2.0*etot);
    x2 = -0.5*ltot*ltot/etot;
    
    if (etot<0.0) {
        rmin = x1-sqrt(x1*x1-x2);
    } else {
        rmin = x1+sqrt(x1*x1-x2);
    }
    
    //      printf("i: %i r: %lf minr: %lf\n",i,r,rmin);
    
    if (rmin<rcap) {
        printf("Capture of particle: %7i time: %10.5lf r: %10.6lf  rmin: %le esmbh: %lf\n",i,t0[i],r,rmin,etot);
        fflush(stdout);
        return 1;
    } else
        return 0;
    
}


int check_gc_disruption(int i) {         //mgc & rh in physical units
    
    double r,omega,lmx,lmy,lmz,dphi2;
    double a[3],dilaton,a1,a2,rj,chi,rmin;
    double RI,MG,MGAL;
    a1 = 0.0;
    a2 = 0.0;
    
    r  = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
    
    if (mgc[i] == 0.0 || rh[i] == 0.0) printf("Warning");
    
    RI = rh[i]/scalesize;
    MG = mgc[i]/scalemass;
    
    if (tdis == 1){
        MGAL = 1.0+msmbh0;
        rj = r*pow(MG/(2.0*MGAL),0.3333333);     //Punktmasse
        chi = RI/rj;
        rmin = 1.0*RI;                          // only when SMBH option is used
    } else if (tdis == 2) {
        lmx     = (y[i]*vz[i]-z[i]*vy[i]);
        lmy     = (z[i]*vx[i]-x[i]*vz[i]);
        lmz     = (x[i]*vy[i]-y[i]*vx[i]);
        omega   = sqrt(lmx*lmx+lmy*lmy+lmz*lmz)/(r*r);
        if (r >= 2.0*RI) {
            dilaton = 1.0+RI/r;
            getforce_scf(i, dilaton, a);
            a1      = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            dilaton = 1.0-RI/r;
            getforce_scf(i, dilaton, a);
            a2      = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        }
        else if (r < 2.0*RI) {
            dilaton = 3.0*RI/r;
            getforce_scf(i, dilaton, a);
            a1      = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            dilaton = 1.0*RI/r;
            getforce_scf(i, dilaton, a);
            a2      = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        }
        
        dphi2   = (a1-a2)/(2.0*RI);
        rj      = pow(MG/(omega*omega-dphi2),0.3333333);
        dilaton = 0.0;
        if ((omega*omega-dphi2) < 0.0) rj = -1.0; //no tidal disruption
        chi     = RI/rj;
        rmin    = 1.0*RI;                                                                                  
    } else {
        rmin = 0.0;
        chi = 0.0;
    }
    
    if (chi>xcrit || (r < rmin && msmbh>0.0)) {                                       
        if (boxorb[i] == 1) {                                                           
            printf("Disruption of GC on Box orbit: %7i Disruptions: %i time: %10.5lf physical time [Myr]  %10.5lf r: %10.6lf\n",i,ngcdis,t0[i],t0[i]*tscale,r);                                      
        } else {                                                                          
            printf("Disruption of GC on non Box orbit: %7i Disruptions: %i time: %10.5lf physical time [Myr] %10.5lf  r: %10.6lf\n",i,ngcdis,t0[i],t0[i]*tscale,r);                                      
        }
        
        fflush(stdout);                                                                
        
        return 1;                                                                      
        
    } else {
        return 0;                                                                      
    }
}


void com_correction() {
    
    int j;
    double cmr[6];
    
    printf("\nApplying centre-of-mass correction.\n");
    
    for (j=0; j<6; j++) cmr[j] = 0.0;
    
    for (j=0; j<nbody; j++) {
        cmr[0] += x[j]*mass[j]; 
        cmr[1] += y[j]*mass[j]; 
        cmr[2] += z[j]*mass[j]; 
        cmr[3] += vx[j]*mass[j]; 
        cmr[4] += vy[j]*mass[j]; 
        cmr[5] += vz[j]*mass[j]; 
    } 
    
    for (j=0; j<nbody; j++) {
        x[j] -= cmr[0];
        y[j] -= cmr[1];
        z[j] -= cmr[2];
        vx[j] -= cmr[3];
        vy[j] -= cmr[4];
        vz[j] -= cmr[5];
    }
    
    return;
}

//********************TTL Integrator*******************************

void get_ttl(int i, double xttl[3], double vttl[3]) {
    
    double v1[3], va[3], vn[3], a[3], apn[3], w, omega, eini;
    double deltae,diff,hini;
    double get_potttl();
    void get_forcettl(),get_force_pnttl();
    double xi[3],vi[3];
    int k,rep,inacc=0;
    double efin,xn,yn,zn,vxn,vyn,vzn,tini,wini,r;
    
    xi[0] = xttl[0];
    xi[1] = xttl[1];
    xi[2] = xttl[2];
    vi[0] = vttl[0];
    vi[1] = vttl[1];
    vi[2] = vttl[2];
    
    
    deltae = 0.0;
    hini = H;
    
    wini = -get_potttl(xi);
    eini = -wini+0.5*(vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
    xn=xi[0];        yn=xi[1];        zn=xi[2];
    vxn=vi[0];       vyn=vi[1];       vzn=vi[2];
    tini = t0[i];
    do {
        xi[0] = xn;        xi[1] = yn;        xi[2] = zn; //x at tini
        vi[0] = vxn;       vi[1] = vyn;       vi[2] = vzn; //v at tini
        t0[i] = tini; //t at tini
        w = wini;
        for (k=0;k<3;k++)
            xi[k] += hini/2.0*vi[k]/w; //xi at tini+hini/2
        t0[i] += hini/2.0/w;  //t at tini+hini/2
        omega = -get_potttl(xi);
        get_forcettl(xi, a);
        get_force_pnttl(i, xi, vi, apn);
        for (k=0;k<3;k++)
            vn[k] = vi[k] + hini*(a[k]+apn[k])/omega;//vn at tini+hini
        do {
            for (k=0;k<3;k++) v1[k] = vn[k];
            for (k=0;k<3;k++) va[k] = (vn[k]+vi[k])/2.0;
            get_force_pnttl(i, xi, va, apn);
            for (k=0;k<3;k++)
                vn[k] = vi[k] + hini*(a[k]+apn[k])/omega;
            diff = (v1[0]-vn[0])*(v1[0]-vn[0])+(v1[1]-vn[1])*(v1[1]-vn[1])+(v1[2]-vn[2])*(v1[2]-vn[2]);
        } while (diff>1.E-20);
        deltae= hini*(va[0]*apn[0]+va[1]*apn[1]+va[2]*apn[2])/omega;
        w += hini/2.0/omega*((vn[0]+vi[0])*a[0]+(vn[1]+vi[1])*a[1]+(vn[2]+vi[2])*a[2]);
        vi[0] = vn[0];  vi[1] = vn[1];  vi[2] = vn[2];//vi at tini+hini
        for (k=0;k<3;k++)
            xi[k] += hini/2.0*vi[k]/w; //xi at tini+hini
        t0[i] += hini/2.0/w;
        efin = get_potttl(xi)+0.5*(vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
        if ((fabs(eini+deltae-efin)>ACCURACY) && (hini>H/256.0)) {
            hini = hini/2.0;
            rep = 1;
        } else
            rep = 0;
    } while (rep==1);
    xttl[0] = xi[0];
    xttl[1] = xi[1];
    xttl[2] = xi[2];
    vttl[0] = vi[0];
    vttl[1] = vi[1];
    vttl[2] = vi[2];
}

void  get_forcettl(double _x[], double a[]) {                                      
    
    int n,l,m;
    double phibar[LMAX][NMAX],dphibar[LMAX][NMAX],p[LMAX][LMAX],dp[LMAX][LMAX],c[LMAX][NMAX],cmod[LMAX][NMAX];
    double elm[LMAX][LMAX],flm[LMAX][LMAX],clm[LMAX][LMAX],dlm[LMAX][LMAX];
    double r,xi,cosphi,phi,costheta,theta;
    double ar,atheta,aphi;
    double fac1,fac2,fac3,fac4;
    double xp[3],xgc[3];                                                      
    
    xgc[0] = _x[0];                                                                          
    xgc[1] = _x[1];                                                                          
    xgc[2] = _x[2];                                                                          
    
    
    r    = sqrt(xgc[0]*xgc[0]+xgc[1]*xgc[1]+xgc[2]*xgc[2]);
    xi   = (r-1.0)/(r+1.0);
    cosphi = xgc[0]/sqrt(xgc[0]*xgc[0]+xgc[1]*xgc[1]);
    phi = acos(cosphi);
    if (xgc[1]<0.0) phi=2.0*PI-phi;
    costheta = xgc[2]/r;
    theta = acos(costheta);
    
    for (l=0;l<LMAX;l++) {
        fac1 = 2.0*(2.0*l+1.5);
        c[l][0] = (double) 1.0;
        c[l][1] = (double) fac1*xi;
        for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*((2.0*n-2.0+fac1)*xi*c[l][n-1]-(n+fac1-2.0)*c[l][n-2]);
    }
    
    for (l=0;l<LMAX;l++) {
        fac1 = (double) -sqrt(4.0*PI)*pow(r,l)/pow(1.0+r,2.0*l+1.0);
        for (n=0;n<NMAX;n++) phibar[l][n]=fac1*c[l][n]*abar[l][n];
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            elm[l][m] = 0.0;
            flm[l][m] = 0.0;
            clm[l][m] = 0.0;
            dlm[l][m] = 0.0;
        }
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) { 
            fac3 = 0.0;
            fac4 = 0.0;
            fac2 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac1 = phibar[l][n];
                fac3 += fac1*clmfull[l][m][n];
                fac4 += fac1*slmfull[l][m][n];
            }
            clm[l][m] += fac3*fac2; 
            dlm[l][m] += fac4*fac2; 
        }
    }
    
    for (l=0;l<LMAX;l++) {
        cmod[l][0] = 0.0;
        fac1 = 1.0/((4.0*l+3.0)*(1.0-xi*xi));
        fac2 = 4.0*l+2.0;
        for (n=1;n<NMAX;n++) cmod[l][n]= (double) fac1*((fac2+n)*c[l][n-1]/c[l][n]-n*xi);   // from Wolfram Math World
    }
    
    
    
    for (l=0;l<LMAX;l++) {
        fac1 = (double) l/r-(2.0*l+1.0)/(1.0+r);
        fac2 = 4.0*(2.0*l+1.5)/pow(1.0+r,2.0);
        for (n=0;n<NMAX;n++) dphibar[l][n] = phibar[l][n]*(fac1+fac2*cmod[l][n]);
    }
    
    
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac3 = 0.0;
            fac4 = 0.0;
            fac2 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac1 = dphibar[l][n];
                fac3 += fac1*clmfull[l][m][n];
                fac4 += fac1*slmfull[l][m][n];
            }
            elm[l][m] += fac3*fac2; 
            flm[l][m] += fac4*fac2; 
        }
    }
    
    fac1 = sqrt(1.0-costheta*costheta);
    
    for (m=0;m<LMAX;m++) {
        p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(fac1,m);
        if (m<LMAX-1) p[m+1][m] = (double) costheta*(2.0*m+1.0)*p[m][m];
        for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costheta*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
    }
    
    fac1 = costheta/sin(theta);
    
    for (l=0;l<LMAX;l++) {
        dp[l][l] = (double) l*fac1*p[l][l]; 
        for (m=0;m<l;m++) dp[l][m] = (double) m*fac1*p[l][m]+p[l][m+1]; 
    }
    
    ar = 0.0; aphi = 0.0; atheta = 0.0;
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac1 = cos(m*phi);
            fac2 = sin(m*phi); 
            ar     += -p[l][m]*(elm[l][m]*fac1+flm[l][m]*fac2);
            atheta += -dp[l][m]/r*(clm[l][m]*fac1+dlm[l][m]*fac2);
            aphi   += -m*p[l][m]/r/sin(theta)*(dlm[l][m]*fac1-clm[l][m]*fac2);
        }
    }
    
    a[0] = sin(theta)*cos(phi)*ar+cos(theta)*cos(phi)*atheta-sin(phi)*aphi;
    a[1] = sin(theta)*sin(phi)*ar+cos(theta)*sin(phi)*atheta+cos(phi)*aphi;
    a[2] = cos(theta)*ar-sin(theta)*atheta;
    
    // ADD SMBH contribution
    xp[0] = xgc[0];  xp[1] = xgc[1]; xp[2] = xgc[2];
    get_force_smbh(xp, a); 
    
    
}   

void get_force_pnttl(int i, double xdf[], double vdf[], double a[]) {
    
    int k;
    double dis,r, vtilde[3], costheta[3];
    double adf[3] = {0.0,0.0,0.0};
    double gamma[3], gammac, vgc,lamb,sigmai,chiparam,bracketp;
    
    if (mgc[i] > 0.0 && t0[i]*tscale > tinit) {   //Calculate vec(a_DF) only for GCs 
        
        dis = sqrt(xdf[0]*xdf[0]+xdf[1]*xdf[1]+xdf[2]*xdf[2]);
        vgc = sqrt(vdf[0]*vdf[0]+vdf[1]*vdf[1]+vdf[2]*vdf[2]);
        r   = sqrt(xdf[0]*xdf[0]+xdf[1]*xdf[1]+xdf[2]*xdf[2]);
        
        if (rh[i] > (mgc[i]*scalesize/scalemass)) {
            lamb = log(dis*scalesize/rh[i]);
        } else {
            lamb = log(dis*scalemass/mgc[i]);
        }
        if (lamb < 0.0) lamb = 0.0;       
        //lamb=6.0; //Be careful! Only for testing!!!
        double rtemp;
        double sigma_1_inter = 0.0;
        double sigma_2_inter = 0.0;
        double sigma_3_inter = 0.0;
        double e1_inter[3] = {0.0, 0.0, 0.0};
        double e2_inter[3] = {0.0, 0.0, 0.0};
        double e3_inter[3] = {0.0, 0.0, 0.0};
        double weight = 0.0;
        
        for (k=0;k<98*GRIDDEPTH+1;k++) {
            if (grid_cell_occupation[k]>2) {
                //squared radial distance of GC from grid-box-COM
                rtemp = sqrt(pow(xdf[0]-grid_coordinates[k][0],2)+pow(xdf[1]-grid_coordinates[k][1],2)+pow(xdf[2]-grid_coordinates[k][2],2));
                if (rtemp < dis) {
                    rtemp = pow(rtemp,64.0);
                    sigma_1_inter += sigma_1[k]/rtemp;
                    sigma_2_inter += sigma_2[k]/rtemp;
                    sigma_3_inter += sigma_3[k]/rtemp;
                    e1_inter[0] += e1[k][0]/rtemp;
                    e1_inter[1] += e1[k][1]/rtemp;
                    e1_inter[2] += e1[k][2]/rtemp;
                    e2_inter[0] += e2[k][0]/rtemp;
                    e2_inter[1] += e2[k][1]/rtemp;
                    e2_inter[2] += e2[k][2]/rtemp;
                    e3_inter[0] += e3[k][0]/rtemp;
                    e3_inter[1] += e3[k][1]/rtemp;
                    e3_inter[2] += e3[k][2]/rtemp;
                    weight += 1.0/rtemp;
                }
            }
        }
        if (weight>0.0) {
            sigma_1_inter /= 1.0*weight;
            sigma_2_inter /= 1.0*weight;
            sigma_3_inter /= 1.0*weight;
            e1_inter[0] /= 1.0*weight;
            e1_inter[1] /= 1.0*weight;
            e1_inter[2] /= 1.0*weight;
            e2_inter[0] /= 1.0*weight;
            e2_inter[1] /= 1.0*weight;
            e2_inter[2] /= 1.0*weight;
            e3_inter[0] /= 1.0*weight;
            e3_inter[1] /= 1.0*weight;
            e3_inter[2] /= 1.0*weight;
        } else {
            sigma_1_inter = 1.0E10;
            sigma_2_inter = 1.0E10;
            sigma_3_inter = 1.0E10;
            e1_inter[0] = 1.0;
            e1_inter[1] = 0.0;
            e1_inter[2] = 0.0;
            e2_inter[0] = 0.0;
            e2_inter[1] = 1.0;
            e2_inter[2] = 0.0;
            e3_inter[0] = 0.0;
            e3_inter[1] = 0.0;
            e3_inter[2] = 1.0;
        }
        
        
        if (dfmethod == 1) { 
            // standard Chandrasekhar
            
            sigmai = (1.0/sqrt(3.0))*sqrt(sigma_1_inter*sigma_1_inter+
                                          sigma_2_inter*sigma_2_inter+
                                          sigma_3_inter*sigma_3_inter); 
            chiparam = vgc/(1.414213562*sigmai);
            gammac = 12.566370616*lamb*mgc[i]*get_densiti(i)/(scalemass*pow(vgc,3)); 
            bracketp = erf(chiparam)-2.0*chiparam/1.772453851*exp(-(chiparam*chiparam));
            for (k=0;k<3;k++) {
                adf[k] = -gammac*bracketp*vdf[k];
                // adf[k] = -0.0005524446*vdf[k]/(pow(dis,2));//1e7!! //For SIS
            }
            
        } else if (dfmethod == 2) {
            // triaxial generalization of Chandrasekhar
            
            //transform velocity vector to eigensystem of corresponding cell 
            
            costheta[0] = (e1_inter[0]*vdf[0]+e1_inter[1]*vdf[1]+e1_inter[2]*vdf[2])/(sqrt(e1_inter[0]*e1_inter[0]+e1_inter[1]*e1_inter[1]+e1_inter[2]*e1_inter[2])*vgc);
            costheta[1] = (e2_inter[0]*vdf[0]+e2_inter[1]*vdf[1]+e2_inter[2]*vdf[2])/(sqrt(e2_inter[0]*e2_inter[0]+e2_inter[1]*e2_inter[1]+e2_inter[2]*e2_inter[2])*vgc);
            costheta[2] = (e3_inter[0]*vdf[0]+e3_inter[1]*vdf[1]+e3_inter[2]*vdf[2])/(sqrt(e3_inter[0]*e3_inter[0]+e3_inter[1]*e3_inter[1]+e3_inter[2]*e3_inter[2])*vgc);
            
            vtilde[0] = vgc*costheta[0];
            vtilde[1] = vgc*costheta[1];
            vtilde[2] = vgc*costheta[2];
            
            double sigma2[3];
            
            sigma2[0] = sigma_1_inter*sigma_1_inter;
            sigma2[1] = sigma_2_inter*sigma_2_inter;
            sigma2[2] = sigma_3_inter*sigma_3_inter;
            sigmai = (1.0/sqrt(3.0))*sqrt(sigma_1_inter*sigma_1_inter+sigma_2_inter*sigma_2_inter+sigma_3_inter*sigma_3_inter); 
            for (k=0;k<3;k++) {
                
                gammac = 5.0132565493*lamb*mgc[i]*get_densiti(i)/(scalemass*pow(sigma_1_inter,3));// sigma_1_inter,3) original
                
                gamma[k] = gammac*get_integral(i, k, vtilde, sigma2); //Attention: sigma_1>sigma_2>sigma_3 as well as eigenvalues e1>e2>e3
            }
            //printf("\ngammas: %f\t%f\t%f\n", gamma[0], gamma[1], gamma[2]);
            
            for (k=0;k<3;k++) {
                //			adf[k] = -gamma[0]*vtilde[0]*e1[grid_cell_nr[i]][k]-gamma[1]*vtilde[1]*e2[grid_cell_nr[i]][k]-gamma[2]*vtilde[2]*e3[grid_cell_nr[i]][k];
                adf[k] = -gamma[0]*vtilde[0]*e1_inter[k]-gamma[1]*vtilde[1]*e2_inter[k]-gamma[2]*vtilde[2]*e3_inter[k];
            }
            
        } 	
        
   	} //end if
    
    a[0] = adf[0];
    a[1] = adf[1];
    a[2] = adf[2];
    
    
}   




double get_potttl(double _x[]) {
    
    int l,m,n;  
    double clm[LMAX][LMAX],dlm[LMAX][LMAX];
    double r,xi,phi,cosphi,costheta,pot;
    double phibar[LMAX][NMAX],p[LMAX][LMAX],c[LMAX][NMAX];
    double fac1,fac2;
    
    r    = sqrt(_x[0]*_x[0]+_x[1]*_x[1]+_x[2]*_x[2]);
    xi   = (r-1.0)/(r+1.0);
    cosphi = _x[0]/sqrt(_x[0]*_x[0]+_x[1]*_x[1]);
    phi = acos(cosphi);
    if (_x[1]<0.0) phi=2.0*PI-phi;
    costheta = _x[2]/r;
    
    for (l=0;l<LMAX;l++) {
        fac1 = 2.0*l+1.5;
        c[l][0] = (double) 1.0;
        c[l][1] = (double) 2.0*fac1*xi;
        for (n=2;n<NMAX;n++) c[l][n] = (double) 1.0/(1.0*n)*(2.0*(n-1.0+fac1)*xi*c[l][n-1]-(n+2.0*fac1-2.0)*c[l][n-2]);
    }
    
    fac1 = -sqrt(4.0*PI)/(1.0+r);
    
    for (l=0;l<LMAX;l++) {
        fac2 = fac1*pow(r,l)/pow(1.0+r,2.0*l);
        for (n=0;n<NMAX;n++) phibar[l][n]=fac2*c[l][n];
    }
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            clm[l][m] = 0.0;
            dlm[l][m] = 0.0;
        }
    }
    
    
    for (l=0;l<LMAX;l++) {
        for (m=0;m<=l;m++) {
            fac1 = narr[l][m];
            for (n=0;n<NMAX;n++) {
                fac2 = fac1*abar[l][n]*phibar[l][n];
                clm[l][m] += fac2*clmfull[l][m][n];
                dlm[l][m] += fac2*slmfull[l][m][n];
            }
        }
    }
    
    
    for (m=0;m<LMAX;m++) {
        p[m][m] = (double)  pow(-1.0,m)*factrl(2*m)/factrl(m)/pow(2.0,m)*pow(1.0-costheta*costheta,0.5*m);
        if (m<LMAX-1) p[m+1][m] = (double) costheta*(2.0*m+1.0)*p[m][m];
        for (l=2+m;l<LMAX;l++) p[l][m] = (double) 1.0/(l-m)*(costheta*(2.0*l-1.0)*p[l-1][m]-(l+m-1.0)*p[l-2][m]);
    }
    
    pot = 0.0;
    
    for (l=0;l<LMAX;l++)
        for (m=0;m<=l;m++)
            pot += (double) p[l][m]*(clm[l][m]*cos(m*phi)+dlm[l][m]*sin(m*phi));
    
    
    // ADD SMBH contribution
    pot += (-1.0*msmbh)/r;
    //  printf("i: %i n: %lf\n",i,pot);
    //  fflush(stdout);
    
    return(pot);
    
}

