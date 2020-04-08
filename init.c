/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Hydrodynamic jet propagation in 2D cylindrical coordinates.

  This problem considers the propagation of a hydrodynamic jet into a
  static uniform medium with constant density and pressure.
  The ambient density, in units of the jet density, is prescribed to be
  \f$ \rho_a = \eta \f$ where \f$\eta\f$ is the ambient/jet density ratio.
  The ambient pressure is \f$ p = 1/\Gamma \f$ when an \c IDEAL EoS is used
  or it is computed from temperature as \f$ p = p(\rho_a, T_a) \f$ for the
   \c PVTE_LAW EoS (here Ta is the ambient temperature).
  These values are set in Init() while the jet inflow is set through 
  a user-defined boundary condition at the lower z-boundary.
  A simple top-hat injection nozzle is used.
  
  The configuration is defined in terms of the following parameters:

  - <tt>g_inputParam[ETA]</tt>:   density ratio between ambient and jet;
  - <tt>g_inputParam[MACH]</tt>:  jet Mach number;
  - <tt>g_inputParam[TJET]</tt>:  jet temperature (only for \c PVTE_LAW EoS).

  defined in \c pluto.ini.
  The reference density and length are given by the jet density and radius
  while the reference velocity is 1 Km/s.
  The actual numerical values are needed only when using the \c PVTE_LAW EoS.
     
  - Configuration #01 uses an \c IDEAL EoS
  - Configurations #02 and #03 set a highly supersonic molecular jet
    evolving with the \c PVTE_LAW EoS.
    The first one adopts the root-finder version while the second one
    adopts the tabulated version.

  \image html hd_jet.jpg "Pressure (left) and density (right) maps for configuration #01 at t=15" 
  \author A. Mignone (mignone@ph.unito.it)
  \date   April 13, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <math.h>

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double Ta = 50.0;  /* Ambient temperature */

  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  #if EOS == IDEAL
   g_gamma = 5.0/3.0;
   v[PRS]  = 1.0/g_gamma;
  #elif EOS == PVTE_LAW
   v[PRS] = Pressure(v,Ta);
  #endif

}
/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 **************************************************************** */
{ }
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions in the lower boundary ghost
 *  zones.  The profile is top-hat: 
 *  \f[
 *     V_{ij} = \left\{\begin{array}{ll}
 *     V_{\rm jet} & \quad\mathrm{for}\quad r_i < 1 \\ \noalign{\medskip}
 *     \mathrm{Reflect}(V)  & \quad\mathrm{otherwise}
 *    \end{array}\right.
 *  \f]
 * where \f$ V_{\rm jet} = (\rho,v,p)_{\rm jet} = (1,M,1/\Gamma)\f$ and
 * \c M is the flow Mach number (the unit velocity is the jet sound speed, 
 * so \f$ v = M\f$).
 *
 *********************************************************************** */
{
  double xbh1, ybh1, zbh1, xbh2, ybh2, zbh2, mbh1 = 2, mbh2 = 1, rbh1, rbh2, d1, d2;
  double r_0 = 1.0, omega = 1, *x1, *x2, *x3, *dx1, *dx2, *dx3;
  int i,j,k;
  double mgas1, mgas2, mdot1, mdot2, modt, dv;
  FILE *file;
  static double saveTime;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  dx1 = grid[IDIR].dx;
  dx2 = grid[JDIR].dx;
  dx3 = grid[KDIR].dx;

  rbh1 = (r_0*mbh2/(mbh1+mbh2));
  rbh2 = (r_0*mbh1/(mbh1+mbh2));

  xbh1 = rbh1*cos(omega*g_time);
  ybh1 = rbh2*sin(omega*g_time);
  zbh1 = 0;
  
  xbh2 = rbh2*cos(omega*g_time+CONST_PI);
  ybh2 = rbh2*sin(omega*g_time+CONST_PI);
  zbh2 = 0;


  mgas1 = 0;
  mgas2 = 0;

  mdot1 = 0;
  mdot2 = 0;
  mdot  = 0;

  if (side == 0){
    TOT_LOOP(k,j,i){
      d1 = sqrt((x1[i] - xbh1)*(x1[i] - xbh1) + (x2[j] - ybh1)*(x2[j] - ybh1) +  (x3[k] - zbh1)*(x3[k] - zbh1));//distance of each cell from bh1
      d2 = sqrt((x1[i] - xbh2)*(x1[i] - xbh2) + (x2[j] - ybh2)*(x2[j] - ybh2) +  (x3[k] - zbh2)*(x3[k] - zbh2));//distance of each cell from bh2
      
      

      if  (d1<0.3)
	{
	  //resets the values of the grid around our defined boundary
	  d->Vc[PRS][k][j][i] = 1/g_gamma;
	  d->Vc[RHO][k][j][i] = 1.0;
	  // d->Vc[VX1][k][j][i] = 0;
	  // d->Vc[VX2][k][j][i] = 0;
	  // d->Vc[VX3][k][j][i] = 0;
	  // d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
	 
	  //calculates our gas mass at our bounary
	  dv = dx1[i]*dx2[j]*dx3[k];
	  mgas1 = mgas1 +  d->Vc[RHO][k][j][i]*dv;
	  
	  //mdot
	  mdot
	}

      if (d2<0.3)
	{
	  //resets the values of the grid around our defined boundary
	  d->Vc[PRS][k][j][i] = 1/g_gamma;
	  d->Vc[RHO][k][j][i] = 1.0;
	  // d->Vc[VX1][k][j][i] = 0;
	  // d->Vc[VX2][k][j][i] = 0;
	  // d->Vc[VX3][k][j][i] = 0;
	  // d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
	  
	  //calculates our gass mass at our bounary
	  dv = dx1[i]*dx2[j]*dx3[k]; 
	  mgas2 = mgas2 +  d->Vc[RHO][k][j][i]*dv;
	}

}
 }
  if(saveTime != g_time)//stops time doubling in output
    {
      saveTime = g_time;
      file = fopen("mdot6", "a");
      fprintf (file, "%12.5f  %12.5f  %12.5f \n", g_time, mgas1, mgas2);
      fclose(file);
    }
}



#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the graviational potential as function of the coordinates.
 *
 *********************************************************************** */
{
  double r, com;
  double xbh1, ybh1, zbh1, xbh2, ybh2, zbh2, mbh1 = 2, mbh2 = 1, rbh1, rbh2, dbh1, dbh2;
  double r_0 = 1.0, omega = 1, Phi1, Phi2;

  r = sqrt(x1*x1 + x2*x2 + x3*x3);

  rbh1 = (r_0*mbh2/(mbh1+mbh2));
  rbh2 = (r_0*mbh1/(mbh1+mbh2));

  xbh1 = rbh1*cos(omega*g_time);
  ybh1 = rbh2*sin(omega*g_time);
  zbh1 = 0;

  dbh1 = (xbh1 - x1)*(xbh1 - x1) + (ybh1 - x2)*(ybh1 - x2) + (zbh1 - x3)*(zbh1 - x3);
  dbh1 = sqrt(dbh1);
  Phi1 = -mbh1/dbh1;

  
  xbh2 = rbh2*cos(omega*g_time+CONST_PI);
  ybh2 = rbh2*sin(omega*g_time+CONST_PI);
  zbh2 = 0;

  dbh2 = (xbh2 - x1)*(xbh2 - x1) + (ybh2 - x2)*(ybh2 - x2) + (zbh2 - x3)*(zbh2 - x3);
  dbh2 = sqrt(dbh2);
  Phi2 = -mbh2/dbh2;

  //com = (mbh1*(xbh1+ybh1) + mbh2*(xbh2+ybh2)) / mbh1+mbh2;

  //xbh2 = r_0*cos(omega*g_time + CONST_PI);
  //ybh2 = r_0*sin(omega*g_time + CONST_PI);
  //zbh2 = 0;

  return Phi1+Phi2;
}
#endif


