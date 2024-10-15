/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
scalar f[];

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  /*
  Actual run and codes!
  */
  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});

  // calculate the center of mass of f[]

  double xcom = 0.0, ucom = 0.0, wt = 0.0;

  foreach() {
    xcom += (2*pi*y*clamp(f[], 0.0, 1.0))*x*sq(Delta);
    ucom += (2*pi*y*clamp(f[], 0.0, 1.0))*u.x[]*sq(Delta);
    wt += (2*pi*y*clamp(f[], 0.0, 1.0))*sq(Delta);
  }

  if (wt > 0.0) {
    xcom /= wt;
    ucom /= wt;
  }

  fprintf(ferr, "%g %g %g", t, xcom, ucom);

}
