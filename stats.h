/** @file stats.h
*   @brief Contains functions for statistical analysis.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date October 2015
*  @bug No known bugs
*/

#ifndef STATS_H
#define STATS_H

//================================================================================================//
//===================================STANDARD INCLUDES============================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"
#include "random_numbers.h"

//================================================================================================//
//====================================DATA STRUCTURES=============================================//
//================================================================================================//

//================================================================================================//
/** @struct histogram_1d_t
*   @brief This structure is the typedef for the histogram_1d_t object.
*/
//================================================================================================//
typedef struct histogram_1d_s histogram_1d_t;
typedef struct histogram_1d_s{
	unsigned int num_bins;
	unsigned int counts[MAX_NUM_BINS];
	double start_bin;
	double end_bin;
	double bin_width;
	double cdf[MAX_NUM_BINS];
	double probabilities[MAX_NUM_BINS];
	double bins[MAX_NUM_BINS];
} histogram_1d_t;


//================================================================================================//
/** @struct histogram_2d_t
*   @brief This structure is the typedef for the histogram_2d_t object.
*/
//================================================================================================//
typedef struct histogram_2d_s histogram_2d_t;
typedef struct histogram_2d_s{
	unsigned int num_x_bins;
	unsigned int num_y_bins;
	unsigned int counts[MAX_NUM_BINS*MAX_NUM_BINS];
	double x_bin_width;
	double y_bin_width;
	double x_bin_range[2];
	double y_bin_range[2];
	double probabilities[MAX_NUM_BINS*MAX_NUM_BINS];
	double x_bins[MAX_NUM_BINS];
	double y_bins[MAX_NUM_BINS];
} histogram_2d_t;


//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//


//================================================================================================//
/**
* @brief This function calculates the gaussian function for double precision data.
*
* @param[in] double x
* @param[in] double mean
* @param[in] double variance
*
* @return double gaussian(x)
*/
//================================================================================================//
double compute_gaussian_dbl(double, double, double);


//================================================================================================//
/**
* @brief This function calculates the two sided gamma function for double precision data.
*
* @param[in] double x
* @param[in] double sigma
*
* @return double two_sided_gamma(x)
*/
//================================================================================================//
double compute_two_sided_gamma_dbl(double,double);


//================================================================================================//
/**
* @brief This function computes the empirical cdf of double precision data.
*
* @param[in] double* values
* @param[in] unsigned int length
* @param[out] double* cdf
*
* @return NONE
*/
//================================================================================================//
void compute_empirical_cdf_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes a gaussian cdf given an array of support.
*
* @param[in] double mean
* @param[in] double variance
* @param[in] unsigned int length
* @param[in] double* support
* @param[out] double* cdf
*
* @return NONE
*/
//================================================================================================//
void compute_gaussian_cdf_dbl(double,double,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function computes a laplace cdf given an array of support.
*
* @param[in] double mu
* @param[in] double b
* @param[in] unsigned int length
* @param[in] double* support
* @param[out] double* cdf
*
* @return NONE
*/
//================================================================================================//
void compute_laplacian_cdf_dbl(double,double,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function computes a two sided gamma cdf given an array of support.
*
* @param[in] double sigma
* @param[in] unsigned int length
* @param[in] double* support
* @param[out] double* cdf
*
* @return NONE
*/
//================================================================================================//
void compute_two_sided_gamma_cdf_dbl(double,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function computes the Kolomogorob-Smirnov Goodness Of Fit Test Statistic.
*
* @param[in] double* empirical_cdf
* @param[in] double* true_cdf
* @param[in] unsigned int length
*
* @return double ks_stat
*/
//================================================================================================//
double compute_ks_statistic_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function creates a 2D histogram.
*
* @param[in] unsigned int num_x_bins
* @param[in] unsigned int num_y_bins
* @param[in] double* x_range
* @param[in] double* y_range
*
* @return histogram_2d_t* self
*/
//================================================================================================//
histogram_2d_t* create_histogram_2d(unsigned int, unsigned int, double*, double*);
						  

//================================================================================================//
/**
* @brief This function adds data to a 2D histogram.
*
* @param[in,out] histogram_2d_t* self
* @param[in] double complex* points
* @param[in] unsigned int num_points
*
* @return NONE
*/
//================================================================================================//
void add_to_histogram_2d(histogram_2d_t*, double complex*, unsigned int);


//================================================================================================//
/**
* @brief This function makes the histogram into a binary data field based on a threshold of counts.
*
* @param[in,out] histogram_2d_t* self
* @param[in] unsigned int threshold
*
* @return NONE
*/
//================================================================================================//
void make_histogram_binary_2d(histogram_2d_t*, unsigned int);


//================================================================================================//
/**
* @brief This function initialized a 1D histogram.
*
* @param[in,out] histogram_1d_t* self
* @param[in] unsigned int num bins
* @param[in] double* bin_range
*
* @return NONE
*/
//================================================================================================//
void initialize_histogram_1d(histogram_1d_t*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a 1D histogram with passed in data.
*
* @param[in,out] histogram_1d_t* self
* @param[in] unsigned int num_bins
* @param[in] double* data
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_histogram_1d_data(histogram_1d_t*, unsigned int, double*, unsigned int);


//================================================================================================//
/**
* @brief This function initializes a kernel density estimation of the passed in data.
*
* @param[in,out] histogram_1d_t* self
* @param[in] unsigned int num_bins
* @param[in] double* data
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_kernel_density_histogram_1d_data(histogram_1d_t*, unsigned int, double*, unsigned int);


//================================================================================================//
/**
* @brief This function adds data to a 1D histogram.
*
* @param[in,out] histogram_1d_t* self
* @param[in] double* data
* @param[in] unsigned int data_length
*
* @return NONE
*/
//================================================================================================//
void add_to_histogram_1d(histogram_1d_t*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function runs a kde on the data and adds the result to the histogram.
*
* @param[in,out] histogram_1d_t* self
* @param[in] double* data
* @param[in] unsigned int data_length
*
* @return NONE
*/
//================================================================================================//
void add_to_histogram_1d_kde(histogram_1d_t*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes the probability array by copying and normalizing the count.
*
* @param[in,out] histogram_1d_t* self
*
* @return NONE
*/
//================================================================================================//
void initialize_histogram_probabilities(histogram_1d_t*);


//================================================================================================//
/**
* @brief This function prints the bin centers, probabilities, and counts.
*
* @param[in] histogram_1d_t* self
* @param[in] char* filename;
*
* @return NONE
*/
//================================================================================================//
void print_histogram_1d(histogram_1d_t*,char*);


//================================================================================================//
/**
* @brief This function computes the entropy of the random variable described by the 1D histogram.
*
* @param[in] histogram_1d_t* self
*
* @return double entropy
*/
//================================================================================================//
double compute_histogram_1d_entropy(histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the sample kullback leibler divergence for two histograms.
*
* @param[in] histogram_1d_t* true_histogram
* @param[in] histogram_1d_t* approx_histogram
*
* @return double divergence
*/
//================================================================================================//
double compute_sample_kullback_leibler_divergence_1d(histogram_1d_t*, histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the smoothed jensen shannon divergence for two histograms.
*
* @param[in] histogram_1d_t* true_histogram
* @param[in] histogram_1d_t* approx_histogram
*
* @return double divergence
*/
//================================================================================================//
double compute_jensen_shannon_divergence_1d(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the angle between two vectors.
*
* @param[in] double* signal1
* @param[in] double* signal2
* @param[in] unsigned int length
*
* @return double correlation_coefficient
*/
//================================================================================================//
double compute_correlation_coefficient_dbl(double*, double*, unsigned int);


//================================================================================================//
/**
* @brief This function computes size of the intersection minus the union
*
* @param[in] double* signal1
* @param[in] double* signal2
* @param[in] unsigned int length
*
* @return double jaccard_index
*/
//================================================================================================//
double compute_jaccard_index_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes size of the intersection over the size of both sets
*
* @param[in] double* signal1
* @param[in] double* signal2
* @param[in] unsigned int length
*
* @return double sorenson_index
*/
//================================================================================================//
double compute_sorenson_index_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function tests if given sample is draw from a specific distribution.
*
* NOTE: This test puts more weight on the tails of a distribution.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double anderson_darling_statistic
*/
//================================================================================================//
double compute_anderson_darling_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Kuiper distance between two cdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double kuiper_statistic
*/
//================================================================================================//
double compute_kuiper_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Von Mises distance between two cdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double von_mises_statistic
*/
//================================================================================================//
double compute_von_mises_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Hellinger distance between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double hellinger_statistic
*/
//================================================================================================//
double compute_hellinger_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Total Variation distance between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double hellinger_statistic
*/
//================================================================================================//
double compute_total_variation_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Separation distance between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double separation_statistic
*/
//================================================================================================//
double compute_separation_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Chi Squared distance between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double chi_squared_statistic
*/
//================================================================================================//
double compute_chi_squared_statistic_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Bhattacharyya distance between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double bhattacharyya_distance
*/
//================================================================================================//
double compute_bhattacharyya_distance_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Bhattacharyya angle between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double bhattacharyya_angle
*/
//================================================================================================//
double compute_bhattacharyya_angle_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the Kantorovich_metric between two pdfs.
*
* @param[in] histogram_1d_t* empirical_histogram
* @param[in] histogram_1d_t* true_histogram
*
* @return double kantorovich_metric
*/
//================================================================================================//
double compute_kantorovich_metric_dbl(histogram_1d_t*,histogram_1d_t*);


//================================================================================================//
/**
* @brief This function computes the cdf of the histogram values
*
* @param[in] histogram_1d_t* self
*
* @return NONE
*/
//================================================================================================//
void compute_histogram_cdf(histogram_1d_t*);


#endif //STATS_H//
