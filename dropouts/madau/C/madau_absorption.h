/*! \file madau_absorption.h
 *  \brief Function declaration for IGM absorption following Madau 1995.
 *
 *  See Madau 1995, ApJ 441, 18
 *  http://adsabs.harvard.edu/abs/1995ApJ...441...18M
 */
#ifndef MADAU_ABSORPTION_H
#define MADAU_ABSORPTION_H
/*! \fn double *madau_absorption(double *lambda, int n_lambda, double z);
 *  \brief Wavelength-dependent absorption owing to neutral H in the IGM.
 */
double *madau_absorption(double *lambda, int n_lambda, double z);
#endif //MADAU_ABSORPTION_H
