/* Fitting Formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies. */
/* Daniel J. Eisenstein & Wayne Hu, Institute for Advanced Study */

/* There are two primary routines here, one to set the cosmology, the
other to construct the transfer function for a single wavenumber k. 
You should call the former once (per cosmology) and the latter as 
many times as you want. */

/* TFmdm_set_cosm() -- User passes all the cosmological parameters as
arguments; the routine sets up all of the scalar quantites needed 
computation of the fitting formula.  The input parameters are: 

1) omega_matter -- Density of CDM, baryons, and massive neutrinos, in units of the critical density. 
2) omega_baryon -- Density of baryons, in units of critical. 
3) omega_hdm    -- Density of massive neutrinos, in units of critical 
4) degen_hdm    -- (Int) Number of degenerate massive neutrino species 
5) omega_lambda -- Cosmological constant, in units of critical.  
6) hubble       -- Hubble constant, in units of 100 km/s/Mpc 
7) redshift     -- The redshift at which to evaluate */

/* TFmdm_onek_mpc() -- User passes a single wavenumber, in units of Mpc^-1.
Routine returns the transfer function from the Eisenstein & Hu
fitting formula, based on the cosmology currently held in the 
internal variables.  The routine returns T_cb (the CDM+Baryon
density-weighted transfer function), although T_cbn (the CDM+
Baryon+Neutrino density-weighted transfer function) is stored
in the global variable tf_cbnu. 
*/

/* We also supply TFmdm_onek_hmpc(), which is identical to the previous
routine, but takes the wavenumber in units of h Mpc^-1. */

/* We hold the internal scalar quantities in global variables, so that
the user may access them in an external program, via "extern" declarations. */

/* Please note that all internal length scales are in Mpc, not h^-1 Mpc! */

float  TFmdm_onek_mpc(float kk);
float TFmdm_onek_hmpc(float kk);

int    TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm, int degen_hdm, float omega_lambda, float hubble, float redshift);

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Convenience from Numerical Recipes in C, 2nd edition */
static float sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/* ------------------------- Global Variables ------------------------ */

/* The following are set in TFmdm_set_cosm() */
float   alpha_gamma,/* sqrt(alpha_nu) */
alpha_nu,/* The small-scale suppression */
beta_c,/* The correction to the log in the small-scale */
num_degen_hdm,/* Number of degenerate massive neutrino species */
f_baryon,/* Baryon fraction */
f_bnu,/* Baryon + Massive Neutrino fraction */
f_cb,/* Baryon + CDM fraction */
f_cdm,/* CDM fraction */
f_hdm,/* Massive Neutrino fraction */
growth_k0,/* D_1(z) -- the growth function as k->0 */
growth_to_z0,/* D_1(z)/D_1(0) -- the growth relative to z=0 */
hhubble,/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
k_equality,/* The comoving wave number of the horizon at equality*/
obhh,/* Omega_baryon * hubble^2 */
omega_curv,/* = 1 - omega_matter - omega_lambda */
omega_lambda_z, /* Omega_lambda at the given redshift */
omega_matter_z,/* Omega_matter at the given redshift */
omhh,/* Omega_matter * hubble^2 */
onhh,/* Omega_hdm * hubble^2 */
p_c,/* The correction to the exponent before drag epoch */
p_cb,/* The correction to the exponent after drag epoch */
sound_horizon_fit,  /* The sound horizon at the drag epoch */
theta_cmb,/* The temperature of the CMB, in units of 2.7 K */
y_drag,/* Ratio of z_equality to z_drag */
z_drag,/* Redshift of the drag epoch */
z_equality;/* Redshift of matter-radiation equality */

/* The following are set in TFmdm_onek_mpc() */
float gamma_eff,/* Effective \Gamma */
      growth_cb,/* Growth factor for CDM+Baryon perturbations */
      growth_cbnu,/* Growth factor for CDM+Baryon+Neutrino pert. */
      max_fs_correction,  /* Correction near maximal free streaming */
      qq,/* Wavenumber rescaled by \Gamma */
qq_eff,/* Wavenumber rescaled by effective Gamma */
qq_nu,/* Wavenumber compared to maximal free streaming */
tf_master,/* Master TF */
tf_sup,/* Suppressed TF */
y_freestream; /* The epoch of free-streaming for a given scale */

/* Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as */
float   tf_cb,/* The transfer function for density-weighted CDM + Baryon perturbations. */
      tf_cbnu;/* The transfer function for density-weighted CDM + Baryon + Massive Neutrino perturbations. */

/* By default, these functions return tf_cb */

/* ------------------------- TFmdm_set_cosm() ------------------------ */
int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm, int degen_hdm, float omega_lambda, float hubble, float redshift)
/* This routine takes cosmological parameters and a redshift and sets up
all the internal scalar quantities needed to compute the transfer function. */
/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos, in units of the critical density. */
/*   omega_baryon -- Density of baryons, in units of critical. */
/*   omega_hdm    -- Density of massive neutrinos, in units of critical */
/*   degen_hdm    -- (Int) Number of degenerate massive neutrino species */
/*        omega_lambda -- Cosmological constant */
/*   hubble       -- Hubble constant, in units of 100 km/s/Mpc */
/*        redshift     -- The redshift at which to evaluate */
/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
sets many global variables for use in TFmdm_onek_mpc() */
{
    float z_drag_b1, z_drag_b2, omega_denom;
    int qwarn;
    qwarn = 0;

theta_cmb = 2.7255/2.7;  // SHIPPED with:  theta_cmb = 2.7280/2.7;/* Assuming T_cmb = 2.728 K */

    /* Look for strange input */
    if (omega_baryon<0.0) {
      fprintf(stderr,
	      "TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
      qwarn = 1;
    }

    if (omega_hdm<0.0) {
      fprintf(stderr,
	      "TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n");
      qwarn = 1;
    }
    if (hubble<=0.0) {
      fprintf(stderr,"TFmdm_set_cosm(): Negative Hubble constant illegal.\n");
      exit(1);  /* Can't recover */
    } 

    else if (hubble>2.0) {
      fprintf(stderr,"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
      qwarn = 1;
    }

    if (redshift<=-1.0) {
      fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
      exit(1);
    } 

    else if (redshift>99.0) {
      fprintf(stderr,
	      "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
      qwarn = 1;
    }

    if (degen_hdm<1) degen_hdm=1;
    num_degen_hdm = (float) degen_hdm;
    /* Have to save this for TFmdm_onek_mpc() */
    /* This routine would crash if baryons or neutrinos were zero, 
       so don't allow that */
    if (omega_baryon<=0) omega_baryon=1e-5;
    if (omega_hdm<=0) omega_hdm=1e-5;

    omega_curv = 1.0-omega_matter-omega_lambda;
    omhh = omega_matter*SQR(hubble);
    obhh = omega_baryon*SQR(hubble);
    onhh = omega_hdm*SQR(hubble);

    f_baryon = omega_baryon/omega_matter;
    f_hdm = omega_hdm/omega_matter;
    f_cdm = 1.0-f_baryon-f_hdm;
    f_cb = f_cdm+f_baryon;
    f_bnu = f_baryon+f_hdm;

    /* Compute the equality scale. */
    z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));/* Actually 1+z_eq */
    k_equality = 0.0746*omhh/SQR(theta_cmb);

    /* Compute the drag epoch and sound horizon */
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag    = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*(1.0+z_drag_b1*pow(obhh,z_drag_b2));
    y_drag    = z_equality/(1.0+z_drag);

    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));

    /* Set up for the free-streaming & infall growth function */
    p_c  = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

    omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+
						  omega_matter*(1.0+redshift));
    omega_lambda_z = omega_lambda/omega_denom;
    omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
    growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/
    (pow(omega_matter_z,4.0/7.0)-omega_lambda_z+
	    (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
    growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)
						       -omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
growth_to_z0 = growth_k0/growth_to_z0;
    
    /* Compute small-scale suppression */
    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
pow(1+y_drag,p_cb-p_c)*
(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
(1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*
(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = sqrt(alpha_nu);
    beta_c = 1/(1-0.949*f_bnu);
    /* Done setting scalar variables */
hhubble = hubble;/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
    return qwarn;
}

/* ---------------------------- TFmdm_onek_mpc() ---------------------- */

float TFmdm_onek_mpc(float kk)
/* Given a wavenumber in Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
/* Input: kk -- Wavenumber in Mpc^-1 */
/* Output: The following are set as global variables:
growth_cb -- the transfer function for density-weighted
CDM + Baryon perturbations. 
growth_cbnu -- the transfer function for density-weighted
CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */
{
    float tf_sup_L, tf_sup_C;
    float temp1, temp2;

    qq = kk/omhh*SQR(theta_cmb);

    /* Compute the scale-dependent growth functions */
    y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*
SQR(num_degen_hdm*qq/f_hdm);
    temp1 = pow(growth_k0, 1.0-p_cb);
    temp2 = pow(growth_k0/(1+y_freestream),0.7);
    growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
    growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

    /* Compute the master function */
    gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
				 (1+SQR(SQR(kk*sound_horizon_fit*0.43))));
    qq_eff = qq*omhh/gamma_eff;

    tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
    tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
    tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));

    qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
    max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/
(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
    tf_master = tf_sup*max_fs_correction;

    /* Now compute the CDM+HDM+baryon transfer functions */
    tf_cb = tf_master*growth_cb/growth_k0;
    tf_cbnu = tf_master*growth_cbnu/growth_k0;

    return tf_cb;
}

/* ---------------------------- TFmdm_onek_hmpc() ---------------------- */

float TFmdm_onek_hmpc(float kk)
/* Given a wavenumber in h Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
/* Input: kk -- Wavenumber in h Mpc^-1 */
/* Output: The following are set as global variables:
growth_cb -- the transfer function for density-weighted
CDM + Baryon perturbations. 
growth_cbnu -- the transfer function for density-weighted
CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */
{
    return TFmdm_onek_mpc(kk*hhubble);
}


int main(){
  int        i;
  double k, Tk;

  FILE*        output;
  char  filepath[200];

  /* Set cosmology. 
  1) omega_matter -- Density of CDM, baryons, and massive neutrinos, in units of the critical density.                                                                                                         
  2) omega_baryon -- Density of baryons, in units of critical.                                                                                                                                                      
  3) omega_hdm    -- Density of massive neutrinos, in units of critical                                                                                                                                              
  4) degen_hdm    -- (Int) Number of degenerate massive neutrino species                                                                                                                                             
  5) omega_lambda -- Cosmological constant, in units of critical.                                                                                                                                                 
  6) hubble       -- Hubble constant, in units of 100 km/s/Mpc                                                                                                                                                 
  7) redshift     -- The redshift at which to evaluate 
  */ 

  //  params = {'om_m': 0.3106, 'om_b': 0.04898, 'om_L': 0.6894, 'h_100': 0.6770, 'sig_8': 0.811, 'Tcmb': 2.7255, 'zscatter': 1100, 'ns': 0.96824, 'As': 2.1073e-9}
  TFmdm_set_cosm(0.3106, 0.04898, 0.0, 3, 0.6894, 0.677, 3.0);

  sprintf(filepath, "EisensteinAndHu.dat");

  output = fopen(filepath, "w");

  fprintf(output, "##  ** Eisenstein and Hu: **");

  fprintf(output, "\n##  k [h Mpc^-1] \t T(k) \t T(k) (cdm + b + mnu).");
  
  for(i=0; i<100000; i++){
    k  = 0.0001 * i;
  
    Tk = TFmdm_onek_hmpc(k);  // Given a wavenumber in h Mpc^-1, return the transfer function for the cosmology held in the global variables.

    fprintf(output, "\n%.6le \t %.6le \t %.6le", k, tf_cb, tf_cbnu);
  }

  fclose(output);

  printf("Done.\n\n");
  
  return 0;
}
