/** @file random_numbers.h
*   @brief Contains functions for generating random variables.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date September 2015
*  @bug No known bugs
*/

#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <complex.h>
#include "macros.h"

//================================================================================================//
//=====================================LOCAL INCLUDES=============================================//
//================================================================================================//

#include "helper.h"
#include "stats.h"


//================================================================================================//
//========================================MACROS==================================================//
//================================================================================================//

#define mask_all_but_highest_bit(u) ((u) & 0x80000000U)
#define mask_all_but_lowest_bit(u) ((u) & 0x00000001U)
#define mask_highest_bit(u) ((u) & 0x7FFFFFFFU)
#define mask_lowest_bit(u) ((u) & 0xFFFFFFFEU)
#define move_highest_bit(u,v) (mask_all_but_highest_bit(u)|mask_highest_bit(v))
#define move_lowest_bit(u,v) (mask_all_but_lowest_bit(u)|mask_lowest_bit(v))



//================================================================================================//
//======================================DATA STRUCTURES===========================================//
//================================================================================================//

//================================================================================================//
/** @struct random_number_generator_t
*   @brief This structure is the typedef for the random_number_generator_t object.
*/
//================================================================================================//
typedef struct random_number_generator_s{
	uint32_t state_array[625];
	uint32_t *next_random_value;	
	uint32_t state_array_length;
	uint32_t period_constant;
	uint32_t magic_constant;
	int32_t left;
	uint32_t initialized;
} random_number_generator_t;



//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//================================================================================================//
/**
* @brief This function seeds the random number generator.
*
* @param[in,out] random_number_generator_t* self 
* @param[in] uint32_t seed
*
* @return NONE
*/
//================================================================================================//
void seed_random_number_generator(random_number_generator_t*, uint32_t);


//================================================================================================//
/**
* @brief This function initializes the random number generator.
*
* @param[in,out] random_number_generator_t* self 
*
* @return NONE
*/
//================================================================================================//
void initialize_random_number_generator(random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates a random number from the generator.
*
* @param[in,out] random_number_generator_t* self 
*
* @return NONE
*/
//================================================================================================//
uint32_t generate_random_number(random_number_generator_t*);


//================================================================================================//
/**
* @brief This function samples from a random gaussian distribution.
*
* @param[in] double mean 
* @param[in] double variance
* @param[in,out] random_number_generator_t* rng 
*
* @return double random_number
*/
//================================================================================================//
double get_random_gaussian(double, double, random_number_generator_t*);


//================================================================================================//
/**
* @brief This function samples from a random uniform distribution.
*
* @param[in] uint32_t low
* @param[in] uint32_t high
* @param[in,out] random_number_generator_t* rng 
*
* @return uint32_t random_number
*/
//================================================================================================//
uint32_t get_random_uniform(uint32_t, uint32_t, random_number_generator_t*);


//================================================================================================//
/**
* @brief This function samples from a random uniform distribution and returns a double value.
*
* @param[in] double low
* @param[in] double high
* @param[in,out] random_number_generator_t* rng 
*
* @return double random_number
*/
//================================================================================================//
double get_random_uniform_dbl(double,double,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precsion samples from a random uniform distribution.
*
* @param[in] double low
* @param[in] double high
* @param[in] unsigned int length
* @param[in,out] process
* @param[in,out] random_number_generator_t* rng 
*
* @return double random_number
*/
//================================================================================================//
void generate_uniform_process_dbl(double,double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random gaussian distribution.
*
* @param[in] double mean
* @param[in] double variance
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_gaussian_process_dbl(double, double, unsigned int, double*, random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random pink gaussian distribution.
*
* @param[in] double mean
* @param[in] double variance
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_pink_gaussian_process_dbl(double,double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random brown gaussian distribution.
*
* @param[in] double mean
* @param[in] double variance
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_brown_gaussian_process_dbl(double,double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random laplacian distribution.
*
* @param[in] double mu
* @param[in] double b
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_laplacian_process_dbl(double,double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random exponential distribution.
*
* @param[in] double gamma
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_exponential_process_dbl(double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random gamma distribution.
*
* @param[in] double k
* @param[in] double theta
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_gamma_process_dbl(double,double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a random beta distribution.
*
* @param[in] double a
* @param[in] double b
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_beta_process_dbl(double,double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates double precision samples from a two sided gamma distribution.
*
* @param[in] double sigma
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_two_sided_gamma_process_dbl(double,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates complex samples from a random gaussian distribution.
*
* @param[in] double real_mean
* @param[in] double imag_mean
* @param[in] double variance
* @param[in] unsigned int length
* @param[in] unsigned int normalize
* @param[in,out] double complex* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_gaussian_process_cmplx(double, double, double, unsigned int, unsigned int, double complex*, random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates samples from a random uniform distribution.
*
* @param[in] uint32_t low
* @param[in] uint32_t high
* @param[in] unsigned int length
* @param[in,out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_uniform_process(uint32_t, uint32_t, unsigned int, uint32_t*, random_number_generator_t*);



//================================================================================================//
/**
* @brief This function randomly samples a double precision array without replacement.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] unsigned int resampled_length
* @param[out] double* resampled
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void randomly_sample_array_dbl(double*,unsigned int,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function randomly samples an unsigned int an array without replacement.
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
* @param[in] unsigned int resampled_length
* @param[out] unsigned int* resampled
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void randomly_sample_array_uint(unsigned int*,unsigned int,unsigned int,
								unsigned int*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function generates an AR(p) process.
*
* @param[in] double* coeffs
* @param[in] unsigned int order
* @param[in] unsigned int length
* @param[out] double* process 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void generate_ar_process_dbl(double*,unsigned int,unsigned int,double*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function adds zero mean gaussian noise to a given vector.
*
* @param[in] double* vector
* @param[in] unsigned int length
* @param[in] double variance
* @param[out] double* noisy_vector 
* @param[in,out] random_number_generator_t* rng 
*
* @return NONE
*/
//================================================================================================//
void add_gaussian_noise_dbl(double*, unsigned int, double, double*, random_number_generator_t*);


//================================================================================================//
/**
* @brief This function creates a korobov sequence of lattice sampled doubles from [0,1]
*
* @param[in] unsigned int base
* @param[in] unsigned int num_samples
* @param[out] double* sequence
*
* @return NONE
*/
//================================================================================================//
void make_korobov_sequence(unsigned int, unsigned int, double* sequence);


//================================================================================================//
/**
* @brief This function generates a bandlimited signal with random frequencies and amplitudes.
*
* @param[in] unsigned int num_frequencies
* @param[in] unsigned int signal_length
* @param[in] double sampling_rate
* @param[in] random_number_generator_t* rng
* @param[out] double* wave
*
* @return NONE
*/
//================================================================================================//
void generate_random_wave(unsigned int,unsigned int,double,random_number_generator_t*,double*);


//================================================================================================//
/**
* @brief This function generates a harmonic signal with random frequencies,amplitudes,and phases.
*
* @param[in] double fundamental
* @param[in] double sampling_rate
* @param[in] unsigned int num_harmonics
* @param[in] unsigned int length
* @param[in] random_number_generator_t* rng
* @param[out] double* wave
*
* @return NONE
*/
//================================================================================================//
void generate_random_harmonic_wave(double,double,unsigned int,unsigned int,
								   random_number_generator_t*,double*);


//================================================================================================//
/**
* @brief This function shuffles an unsigned int array into a random order.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int length
* @param[in] random_number_generator_t* rng
*
* @return NONE
*/
//================================================================================================//
void shuffle_array_uint(unsigned int*,unsigned int,random_number_generator_t*);



//================================================================================================//
/**
* @brief This function tests the array shuffle routine.
*
* @return NONE
*/
//================================================================================================//
void test_shuffle_array();


//================================================================================================//
/**
* @brief This function tests the random process generation routines.
*
* @return NONE
*/
//================================================================================================//
void test_random_processes();


//================================================================================================//
/**
* @brief This function writes the various colors of noise to wav files.
*
* @return NONE
*/
//================================================================================================//
void test_random_noise();

#endif //RANDOM_NUMBERS_H//
