#include "random_numbers.h"


//================================================================================================//
//=======================================RANDOM METHODS===========================================//
//================================================================================================//

void seed_random_number_generator( random_number_generator_t* self, 
								   uint32_t seed )
{
	//===Seed RNG===//
    srand(seed); 

	//===Set True===//
	self->initialized = 1;

	return;
}

void initialize_random_number_generator( random_number_generator_t* self )
{
	uint32_t seed;

	//===Set Fields===//
	self->initialized = 0;

	//===Seed Generator===//
	seed = time(NULL);
    seed_random_number_generator(self,seed);
	
	return;		
}

uint32_t generate_random_number( random_number_generator_t* self )
{

	//===Locals===//    
	uint32_t y;

	//===Check if Initialized===//
	if (self->initialized == 0){
		seed_random_number_generator(self, time(NULL));
	}

	//===Generate Random Number===//
	y = rand();

    return y;
}

double get_random_gaussian( double mean, 
						    double variance,
							random_number_generator_t* rng )
{

	uint32_t random;
	double fac, rsq, v1, v2;

	rsq = 2;
	while (rsq >= 1.0 || fabs(rsq) < 10.0*EPS){
		random = generate_random_number(rng);
		v1 = 2.0*(((double)random)/(double)UINT32_MAX)-1.0;
		random = generate_random_number(rng);
		v2 = 2.0*(((double)random)/(double)UINT32_MAX)-1.0;
		rsq = v1*v1 + v2*v2;
	}
	fac = sqrt(-2.0*log(rsq)/rsq);

	if (fac > sqrt(-2.0*log(0.5)/0.5)){
		fac *= v2;
	}
	else if (fac < sqrt(-2.0*log(0.5)/0.5)){
		fac *= v1;
	}
	else{
		fac *= 0.5*(v1+v2);
	}

	fac = sqrt(variance)*fac + mean;

	return fac;
}

uint32_t get_random_uniform( uint32_t low, 
							 uint32_t high,
							 random_number_generator_t* rng )
{
    uint32_t val;
    uint32_t range;
    uint32_t scale;
	uint32_t swap_temp;
	
	if (low > high){
		SWAP(low, high);
		fprintf(stderr, "Error: Low Is Greater Than High! Swapping! In Function -- get_random_uniform!\n");
	}
	if (low == high){
		fprintf(stderr, "Error: Low Is Equal To High! In Function -- get_random_uniform!\n");
		return 0;
	}
	if (high - low > RAND_MAX){
		fprintf(stderr, "Error: Range Is Too Large! In Function -- get_random_uniform!\n");
		return 0;
	}

	range = high-low;
	scale = ULONG_MAX/range;
	val = scale * range + 1;
	while (val >= scale * range){
        val = generate_random_number(rng);
	}
    return val/scale + low;
}

double get_random_uniform_dbl( double low,
							   double high,
							   random_number_generator_t* rng )
{
	uint32_t val;
	double div;
	double swap_temp;

	if (low > high){
		SWAP(low, high);
		fprintf(stderr, "Error: Low Is Greater Than High! Swapping! In Function -- get_random_uniform_dbl!\n");
	}
	if (fabs(low-high) < 10.0*EPS){
		fprintf(stderr, "Error: Low Is Equal To High! In Function -- get_random_uniform_dbl!\n");
		return 0;
	}
	if (high - low > DBL_MAX){
		fprintf(stderr, "Error: Range Is Too Large! In Function -- get_random_uniform_dbl!\n");
		return 0;
	}

	//===Get Uniform Number===//
	val = get_random_uniform(0, RAND_MAX, rng);
	
	//===Convert To Double===//
    div = RAND_MAX / (high-low);

	return  low + val/div;
}

void generate_uniform_process_dbl( double low,
								   double high,
								   unsigned int length,
								   double* process,
								   random_number_generator_t* rng )
{
	unsigned int i;
	for (i=0; i<length; i++){
		process[i] = get_random_uniform_dbl(low, high, rng);
	}
	return;
} 


void generate_gaussian_process_dbl( double mean,
									double variance,
									unsigned int length,
									double* process,
									random_number_generator_t* rng )
{
	unsigned int i;	
	for (i=0; i<length; i++){
		process[i] = get_random_gaussian(mean, variance, rng);
	}

	return;
}

void generate_laplacian_process_dbl( double mu,
								     double b, 
								 	 unsigned int length,
								 	 double* process,
								 	 random_number_generator_t* rng )
{

	unsigned int i;

	//===Initialize Process===//
	generate_uniform_process_dbl(-0.5, 0.5, length, process, rng);

	//===Make Laplacian===//
	for (i=0; i<length; i++){
		process[i] = mu - b * sign_dbl(process[i]) * log(1.0 - 2.0*fabs(process[i]));
	}
	
	return;
}

void generate_exponential_process_dbl( double gamma,
								       unsigned int length,
								       double* process,
								       random_number_generator_t* rng )
{
	unsigned int i;

	//===Make Uniform===//
	generate_uniform_process_dbl(0.0, 1.0, length, process, rng);

	//===Make Exponential===//
	for (i=0; i<length; i++){
		process[i] = -log(process[i])/gamma;
	}

	return;
}

void generate_gamma_process_dbl( double k,
								 double theta,
								 unsigned int length,
								 double* process,
								 random_number_generator_t* rng )
{
	unsigned int i, j;
	double exponential[MAX_SIGNAL_LENGTH];

	for (i=0; i<length; i++){
		generate_exponential_process_dbl(theta, k, exponential, rng);
		for (j=0; j<k; j++){
			process[i] += exponential[j];
		}		
	}

	return;
}

void generate_beta_process_dbl( double a,
								double b,
								unsigned int length,
								double* process,
								random_number_generator_t* rng )
{
	unsigned int i, good;
	double c, alpha, y, u;
	c = 1.0; alpha = 1.0;

	for (i=0; i<length; i++){

		good = 0;
		while (good == 0){
				
			//===Generate Y===//
			generate_uniform_process_dbl(0, 1, 1, &y, rng);			
						
			//===Generate U===//
			generate_uniform_process_dbl(0, 1, 1, &u, rng);
			
			//===Compare===//
			if (u <= (c * pow(y,a-1) * pow((1.0 - y),b-1))/(alpha*1.0)){
				good = 1;
			}
		}
		process[i] = y;
	}
	
	return;
}

void generate_two_sided_gamma_process_dbl( double sigma,
									   	   unsigned int length,
									       double* process,
									   	   random_number_generator_t* rng)
{
	unsigned int i;
	double uniform, y, sigmag, fy, gy, M;

	//===Generate Process Via Rejection Sampling===//
	i = 0;
	sigmag = 10.0 * sigma;
	M = 250.0; 
	while (i < length){

		//===Generate Gaussian===//
		y = get_random_gaussian(0, sigmag, rng);

		//===Compute f(y) and g(y)===//
		fy = compute_two_sided_gamma_dbl(y, sigma);
		gy = compute_gaussian_dbl(y,0,sigmag);

		//===Generate Uniform===//
		uniform = get_random_uniform_dbl(0, 1, rng);

		//===Add To Process===//
		if (uniform < fy/(M*gy)){
			process[i++] = y;
		}
	}

	return;
}
									   


void generate_gaussian_process_cmplx( double real_mean,
									  double imag_mean,
									  double variance,
									  unsigned int length,
									  unsigned int normalize, 
							  		  double complex* process, 
							 		  random_number_generator_t* rng )
{
	unsigned int i;

	//===Sanity Check===//
	if (normalize == 1){
		fprintf(stderr, "Error:: Cannot Normalize! In Function -- generate_gaussian_process_cmplx!\n");
		quit();
	}

	//===Generate Process===//
	for (i=0; i<length; i++){
		process[i] = get_random_gaussian(real_mean, variance, rng) + I*get_random_gaussian(imag_mean, variance, rng);
	}

	return;
}

void generate_uniform_process( uint32_t low,
							   uint32_t high,
							   unsigned int length,
							   uint32_t* process,
							   random_number_generator_t* rng )
{
	unsigned int i;	
	for (i=0; i<length; i++){
		process[i] = get_random_uniform(low, high, rng);
	}

	return;
}

void randomly_sample_array_dbl( double* array,
						 		unsigned int length,
						 		unsigned int resampled_length,
						 		double* resampled,
						 		random_number_generator_t* rng )
{
	unsigned int i, num_indices;
	unsigned int num_so_far, good_sample, sample;
	unsigned int indices[MAX_POLYNOMIAL_ORDER+1];

	//===Sanity Check===//
	if(length > MAX_POLYNOMIAL_ORDER){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- randomly_sample_array_dbl!\n");
		quit();
	}
	if(resampled_length > length){
		fprintf(stderr, "Error:: Resampled length Is Too Large! In Function -- randomly_sample_array_dbl!\n");
		quit();
	}

	//===Check If The Same Length===//
	if (resampled_length == length){
		copy_array_dbl(array, length, resampled);
	}
	else{
		//===Generate Indices===//
		num_so_far = 0; 
		while(num_so_far < resampled_length){
			//===Generate Sample===//
			if (num_so_far == 0){
				indices[num_so_far++] = get_random_uniform(0, length-1, rng);
			}
			else{
				//===Generate Unique Sample===//
				good_sample = 0;
				while (good_sample == 0){
					//===Generate Sample===//
					sample = get_random_uniform(0, length-1, rng);
					//===See If Exists Already===//
					num_indices = UINT_MAX;
					find_indices_where_equal_uint(indices, num_so_far, sample, &num_indices, NULL);
					if (num_indices == 0 && sample != indices[0]){
						good_sample = 1;
					}
				}
				//===Add To List===//
				indices[num_so_far++] = sample;
			}
		}
		//===Copy Over===//
		for (i=0; i<resampled_length; i++){
			resampled[i] = array[indices[i]];
		}
	}

	return;
}

void randomly_sample_array_uint( unsigned int* array,
						 		 unsigned int length,
						 		 unsigned int resampled_length,
						 		 unsigned int* resampled,
						 		 random_number_generator_t* rng )
{
	unsigned int i, num_indices;
	unsigned int num_so_far, good_sample, sample;
	unsigned int indices[MAX_SIGNAL_LENGTH];

	//===Sanity Check===//
	if(length > MAX_SIGNAL_LENGTH){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- randomly_sample_array_uint!\n");
		quit();
	}
	if(resampled_length > length){
		fprintf(stderr, "Error:: Resampled length Is Too Large! In Function -- randomly_sample_array_uint!\n");
		quit();
	}

	//===Check If The Same Length===//
	if (resampled_length == length){
		copy_array_uint(array, length, resampled);
	}
	else{
		//===Generate Indices===//
		num_so_far = 0; 
		while(num_so_far < resampled_length){
			//===Generate Sample===//
			if (num_so_far == 0){
				indices[num_so_far++] = get_random_uniform(0, length-1, rng);
			}
			else{
				//===Generate Unique Sample===//
				good_sample = 0;
				while (good_sample == 0){
					//===Generate Sample===//
					sample = get_random_uniform(0, length-1, rng);
					//===See If Exists Already===//
					num_indices = UINT_MAX;
					find_indices_where_equal_uint(indices, num_so_far, sample, &num_indices, NULL);
					if (num_indices == 0 && sample != indices[0]){
						good_sample = 1;
					}
				}
				//===Add To List===//
				indices[num_so_far++] = sample;
			}
		}
		//===Copy Over===//
		for (i=0; i<resampled_length; i++){
			resampled[i] = array[indices[i]];
		}
	}

	return;
}


void generate_ar_process_dbl( double* coeffs,
						  	  unsigned int order,
						  	  unsigned int length,
						  	  double* process,
							  random_number_generator_t* rng  )
{
	unsigned int i, j;
	double input;

	//===Create Input Signal===//
	input = get_random_gaussian(0.0, 1.0, rng);

	for (i=0; i<length; i++){
		if (i==0){
			process[i] = input;
		}
		else{
			for (j=1; j<=MIN(i,order); j++){
				process[i] += coeffs[j-1] * process[i-j];
			}
			process[i] += get_random_gaussian(0.0, 1.0, rng);
		}
	}

	return;
}


void add_gaussian_noise_dbl( double* vector,
							 unsigned int length,
							 double variance,
							 double* noisy_vector,
							 random_number_generator_t* rng )
{
	unsigned int i;
	for (i=0; i<length; i++){
		noisy_vector[i] = vector[i] + get_random_gaussian(0.0, variance, rng);
	}	

	return;
}



//================================================================================================//
//==================================PSEUDO RANDOM METHODS=========================================//
//================================================================================================//

void make_korobov_sequence( unsigned int base,
						    unsigned int num_samples, 
							double* sequence )
{
	unsigned int i;

	sequence[0] = fmod(base * 1.0/(double)num_samples, 1.0);
	for (i=1; i<num_samples; i++){
		sequence[i] = fmod(sequence[i-1] + 1.0/(double)num_samples, 1.0);
	}	
	for (i=0; i<num_samples; i++){
		if (sequence[i] < 0.5){
			sequence[i] *= 2.0;
		}
		else{
			sequence[i] = 2.0 * (1.0 - sequence[i]);
		}
	}

	return;
}

void shuffle_array_uint( unsigned int* array,
						 unsigned int length,
						 random_number_generator_t* rng )
{

	int i, j;
	unsigned int temp;

	for (i=0; i<((int)length-1); i++){
		j = i + get_random_uniform(0, length-1-i, rng);
		temp = array[j];
		array[j] = array[i];
		array[i] = temp;
	}

	return;
}

//================================================================================================//
//========================================TEST METHODS============================================//
//================================================================================================//

void test_shuffle_array()
{
	unsigned int i, length;
	unsigned int array[100];
	random_number_generator_t* rng;

	//===Make RNG===//
	rng = malloc(sizeof(random_number_generator_t));
	initialize_random_number_generator(rng);

	//===Create Array===//
	length = 15;
	for (i=0; i<length; i++){
		array[i] = i;
	}
	newline();
	fprintf(stdout, "Original: \n");
	print_vector_uint(array, length, stdout);
	newline();

	//===Shuffle Array===//
	shuffle_array_uint(array, length, rng);
	newline();
	fprintf(stdout, "Shuffled: \n");
	print_vector_uint(array, length, stdout);
	newline();

	//===Clean Up===//
	free(rng);

	return;
}

void test_random_processes()
{

	unsigned int length;
	double process[MAX_SIGNAL_LENGTH];
	random_number_generator_t rng;
	FILE* fout;

	//===Initialize RNG===//
	initialize_random_number_generator(&rng);

	//===Make Laplacian===//
	length = 5000;
	generate_laplacian_process_dbl(0.0, 1.0, length, process, &rng);

	//===Print===//
	fout = fopen("laplacian.dat", "w");
	print_vector_dbl(process, length, fout);
	fclose(fout);

	//===Make Exponential===//
	generate_exponential_process_dbl(1.0, length, process, &rng);

	//===Print===//
	fout = fopen("exponential.dat", "w");
	print_vector_dbl(process, length, fout);
	fclose(fout);

	//===Make Gamma===//
	generate_gamma_process_dbl(1.0, 1.0, length, process, &rng);

	//===Print===//
	fout = fopen("gamma.dat", "w");
	print_vector_dbl(process, length, fout);
	fclose(fout);

	//===Make Beta===//
	generate_beta_process_dbl(2.0, 5.0, length, process, &rng);

	//===Print===//
	fout = fopen("beta.dat", "w");
	print_vector_dbl(process, length, fout);
	fclose(fout);

	//===Make Two Sided Gamma===//
	generate_two_sided_gamma_process_dbl(0.05, length, process, &rng);

	//===Print===//
	fout = fopen("two_sided_gamma.dat", "w");
	print_vector_dbl(process, length, fout);
	fclose(fout);

	return;
}

void test_random_noise()
{
	double *process;
	random_number_generator_t rng;

	//===Mallocs===//
	process = malloc(3*MAX_SIGNAL_LENGTH*sizeof(double));	

	//===Initialize RNG===//
	initialize_random_number_generator(&rng);

	//===Make White Noise===//
	initialize_array_dbl(process, MAX_SIGNAL_LENGTH);
	generate_gaussian_process_dbl(0, 1, MAX_SIGNAL_LENGTH, process, &rng);

	//===Clean Up===//
	free(process);

	return;
}
