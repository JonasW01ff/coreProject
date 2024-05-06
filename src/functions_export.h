/* File: Ratings.h
 * Interface to libRatings.A.dylib 1.0.
 *************************************/
#pragma once
 extern "C" {
/* Adds 'rating' to the set.
 *      rating: Each character adds 1 to the numeric rating
 *              Example: "" = 0, "*" = 1, "**" = 2, "wer " = 4.
 */
void addRating(char *rating);
 
/* Returns the number of ratings in the set.
 */
int ratings(void);
 
/* Returns the mean rating of the set.
 */
char *meanRating(void);
 
/* Clears the set.
 */
void clearRatings(void);

/* Calculates the incomplete Beta function
*/
double xIncomepleteBeta(double a, double b, double x);

/* Calculates the incomplete Beta function
*/
double xBinomialCDF(int k, int n, double p);

/* Calculates alphahage
*/
double xAlphaHage(double R, double n);

/* Function to calculate the annuity rate
*/
double xAnnuityRate(double pD0, double pC, double pn);

/* Function to calculate the annuity length
*/
double xAnnuityLength(double pD0, double pC, double pR);


/* Function to calculate the Confidence Interval for the VaR
*/
int xVaRCI(double *losses, int n, double alpha, double prob, double *res);

/* Function to calculate the MLE of the Gamma distribution
*/
int xGammaMLE(double *bins, short bins_size, int *freqs, int freqs_size, double *alpha, double *beta);
 }


