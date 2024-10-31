/* Title: Coalescence of bubbles

# Version 1.0
# Last modified: Oct 15, 2024

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

// 1 is bubble
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase-tag.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "distance.h"
#include "tag.h"

int MAXlevel; // command line input

#define tsnap (1e-2)
#define tsnap2 (1e-4)
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define VelErr (1e-2)                            // error tolerances in velocity
#define TOL (1e-2)                                 // error tolerance in position

// boundary conditions
f[left] = dirichlet(0.0);
u.t[left] = dirichlet(0.0);

// Other command - line inputs.
double tmax, MuRin, OhOut, RhoIn;
double Rr, zWall;
double Ldomain; // Dimension of the bugger drop and the domain
char nameOut[80], dumpFile[80];

int main(int argc, char const *argv[]) {
  if (argc < 7){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
    return 1;
  }

  // Values taken from the terminal
  MuRin = 1e-2;

  OhOut = atof(argv[1]);
  RhoIn = atof(argv[2]);
  Rr = atof(argv[3]);
  MAXlevel = atoi(argv[4]);
  tmax = atof(argv[5]);
  zWall = atof(argv[6]);

  Ldomain = zWall+2.+2.*Rr+4.0;

  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);

  L0=Ldomain;
  origin(-2.0-zWall, 0.0);
  init_grid (1 << (6));

  rho1 = RhoIn; mu1 = MuRin*OhOut;
  rho2 = 1e0; mu2 = OhOut;
  f.sigma = 1.0;
  ftag.sigma = 0.0;

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "dump");

  run();
}

event init(t = 0){
  if (!restore (file = dumpFile)){
    char filename[60];
    sprintf(filename,"InitialConditionRr-%3.2f.dat", Rr);

    char comm[160];
    sprintf (comm, "scp -r DataFiles/%s .", filename);
    system(comm);

    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8}, MAXlevel).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fraction of the domain. */
    fractions (phi, f);
    fractions (phi, ftag);
    foreach(){
    u.x[] = 0.0;
    u.y[] = 0.0;
    // p[] = 2*(1.-f[]);
    }
    dump (file = dumpFile);
    static FILE * fp2;
    fp2 = fopen("InitialConditionStatus.dat","w");
    fprintf(fp2, "Initial condition is written to %s\n", dumpFile);
    fclose(fp2);
    // return 1;
  }
}

event adapt(i++){
  adapt_wavelet ((scalar *){f, u.x, u.y},
     (double[]){fErr, VelErr, VelErr},
      MAXlevel, MAXlevel-6);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event end (t = end) {
  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, Oh2 %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);
}

scalar posEq[], posPoles[];
event logWriting (t = 0; t += tsnap2; t <= tmax+tsnap) {
  foreach(){
    ftag[] = f[];
  }
  scalar d[];
  // tag all liquid parts starts
  double threshold = 1e-4;
  foreach(){
    d[] = (ftag[] > threshold);
  }

  int n = tag (d), size[n];
  for (int i = 0; i < n; i++){
    size[i] = 0;
  }

  foreach_leaf(serial){
    if (d[] > 0){
      size[((int) d[]) - 1]++;
    }
  }

  #if _MPI
  MPI_Allreduce (MPI_IN_PLACE, size, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif

  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++){
    // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize){
      MaxSize = size[i];
      MainPhase = i+1;
    }
  }
  
  foreach(){
    if(d[] != MainPhase){
      ftag[] = 0.0;
    }
  }
  // tag all liquid parts ends


  double ke = 0., wt = 0., xCOM = 0., Vcm = 0.;

  foreach (reduction(+:ke), reduction(+:wt), reduction(+:xCOM), reduction(+:Vcm)){
    ke += 2*pi*y*(0.5*clamp(ftag[], 0., 1.)*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    xCOM += 2*pi*y*x*clamp(ftag[], 0.0, 1.0)*sq(Delta);
    Vcm += 2*pi*y*u.x[]*clamp(ftag[], 0.0, 1.0)*sq(Delta);
    wt += 2*pi*y*clamp(ftag[], 0.0, 1.0)*sq(Delta);
  }
  xCOM /= wt;

  int nEq = 0;
  double Req = 0.0;
  position (ftag, posEq, {0,1});
  foreach(reduction(+:Req), reduction(+:nEq)){
    if (x < xCOM+TOL && x > xCOM-TOL && posEq[] != nodata){
      Req += posEq[];
      nEq++;
    }
  }

  if (nEq != 0){
    Req /= nEq;
  } else {
    Req = 0.0;
  }

  double zNP = 0.0, zSP = 0.0;
  int nNP = 0, nSP = 0;
  position (ftag, posPoles, {1,0});
  foreach(reduction(+:zNP), reduction(+:nNP), reduction(+:zSP), reduction(+:nSP)){
    if (y < TOL && posPoles[] != nodata){
      if (posPoles[] > 0){
        zNP += posPoles[];
        nNP++;
      } else {
        zSP += posPoles[];
        nSP++;
      }
    }
  }

  if (nNP != 0){
    zNP /= nNP;
  } else {
    zNP = 0.0;
  }
  if (nSP != 0){
    zSP /= nSP;
  } else {
    zSP = 0.0;
  }
  
  static FILE * fp;

  if (pid() == 0) {
    if (i == 0) {
      fprintf (ferr, "i dt t ke Xc Vcm Re ZNp ZSp\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);
      fprintf (fp, "i dt t ke Xc Vcm Re ZNp ZSp\n");
      fprintf (fp, "%d %g %g %g %g %g %g %g %g\n", i, dt, t, ke, xCOM, Vcm/wt, Req, zNP, zSP);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g %g %g %g %g\n", i, dt, t, ke, xCOM, Vcm/wt, Req, zNP, zSP);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g %g %g %g %g\n", i, dt, t, ke, xCOM, Vcm/wt, Req, zNP, zSP);
  }

  assert(ke > -1e-10);
  // dump (file = "dump");
  if (ke < 1e-5 && i > 1000){
    if (pid() == 0){
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen ("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
    }
    return 1;
  }
}
