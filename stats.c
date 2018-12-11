#include "stats.h"

//================================================================================================//
//===================================DISTRIBUTION METHODS=========================================//
//================================================================================================//

double compute_gaussian_dbl( double x,
					 		 double mean,
					 		 double variance )
{
	double constant, exponential;
	constant = 1.0/(variance * sqrt(2.0 * M_PI));
	exponential = exp((-pow((x-mean),2.0))/(2.0*pow(variance,2.0)));
	return constant*exponential;
}

double compute_two_sided_gamma_dbl( double x,
								    double sigma )
{
	double temp;	
	temp = sqrt(sqrt(6.0)/(8.0 * M_PI * sigma * fabs(x)));
	temp *= exp((-sqrt(3.0) * fabs(x))/(sqrt(2.0) * sigma));	
	return temp;
}


void compute_gaussian_cdf_dbl( double mean,
						   	   double variance,
						   	   unsigned int length,
						   	   double* support,
						   	   double* cdf )
{
	unsigned int i;
	double sorted_support[2*(MAX_POLYNOMIAL_ORDER+1)];
	if (length > 2*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- compute_gaussian_cdf_dbl!\n");
		quit();
	}
	
	copy_array_dbl(support, length, sorted_support);
	sort_array_dbl(sorted_support, length);
	for (i=0; i<length; i++){
		cdf[i] = 0.5 * (1.0 + erf( (sorted_support[i] - mean)/(sqrt(variance) * sqrt(2.0))) ); 
	}

	return;
}

void compute_laplacian_cdf_dbl( double mu,
								double b,
								unsigned int length,
								double* support,
								double* cdf )
{

	unsigned int i;
	double sorted_support[2*(MAX_POLYNOMIAL_ORDER+1)];
	if (length > 2*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- compute_gaussian_cdf_dbl!\n");
		quit();
	}
	
	copy_array_dbl(support, length, sorted_support);
	sort_array_dbl(sorted_support, length);
	for (i=0; i<length; i++){

		if (sorted_support[i] < mu){
			cdf[i] = 0.5 * exp((sorted_support[i]-mu)/b);
		}
		else{
			cdf[i] = (1.0 - (0.5 * exp(-1.0 * (sorted_support[i]-mu)/b)));
		}
	}
	return;
}	

void compute_two_sided_gamma_cdf_dbl( double sigma,
									  unsigned int length,
									  double* support,
									  double* cdf )
{
	unsigned int i;
	double x, temp;
	double sorted_support[2*(MAX_POLYNOMIAL_ORDER+1)];
	if (length > 2*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- compute_gaussian_cdf_dbl!\n");
		quit();
	}
	
	copy_array_dbl(support, length, sorted_support);
	sort_array_dbl(sorted_support, length);
	for (i=0; i<length; i++){
		x = sorted_support[i];
		temp = erf(0.5 * ((pow(2.0, 3.0/4.0) * pow(3.0, 1.0/4.0) * sqrt(fabs(x)))/sqrt(sigma)));
		cdf[i] = 0.5 * sqrt(1.0/(sigma*fabs(x))) * sqrt(fabs(x)) * sqrt(sigma) * temp;
		if (x < 0 || fabs(x) < 2.0*EPS){
			cdf[i] *= -1.0;
		}
		cdf[i] += 0.5;
	}
	return;
}


void compute_empirical_cdf_dbl( double* values,
								unsigned int length,
								double* cdf )
{
	unsigned int i, j;
	double sorted[2*(MAX_POLYNOMIAL_ORDER+1)];
	double num_less[2*(MAX_POLYNOMIAL_ORDER+1)];
	if (length > 2*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- compute_empirical_cdf_dbl!\n");
		quit();
	}

	//===Sort Values===//
	copy_array_dbl(values, length, sorted);
	sort_array_dbl(sorted, length);

	//===Find Num Less===//
	for (j=0; j<length; j++){
		num_less[j] = 0;
		for (i=0; i<length; i++){
			if ( sorted[i] < sorted[j] ){
				num_less[j] += 1;
			}
		}
		cdf[j] = num_less[j]/((double)length);
	}	

	return;
}

double compute_ks_statistic_dbl( double* empirical_cdf,
							 	 double* true_cdf,
							 	 unsigned int length )
{
	unsigned int i;
	double diff, maxx;

	maxx = -(DBL_MAX-1.0);
	for (i=0; i<length; i++){
		diff = fabs(empirical_cdf[i] - true_cdf[i]);
		if (diff > maxx){
			maxx = diff;
		}
	}
	return maxx;
}

/*
	When adding things to these histograms, you should only increment counts
	Probability array should be created by copying the counts and then sum normalizing

	Optimal Bin Width:
		optimal_width = 3.0 * 49.0 * var(data)**0.5 * len(data)**(1.0/3.0)
		so then num_bins = (max(data) - min(data))/optimal_width
	From: On Optimal and Data-Based Histograms by David W Scott 1979
*/

//================================================================================================//
//=================================HISTOGRAM 2D METHODS===========================================//
//================================================================================================//

histogram_2d_t* create_histogram_2d( unsigned int num_x_bins,
						  			 unsigned int num_y_bins,
						  			 double* x_range,
						  			 double* y_range )
{
	unsigned int i;
	double x_step_size, y_step_size;
	histogram_2d_t* histogram;
	histogram = malloc(sizeof(histogram_2d_t));

	//===Set Bin Widths===//
	if (num_x_bins % 2 == 0){
		num_x_bins++;		
	}
	if (num_y_bins % 2 == 0){
		num_y_bins++;		
	}
	x_step_size = (x_range[1] - x_range[0])/(double)(num_x_bins-1);	
	histogram->x_bin_width = x_step_size; histogram->num_x_bins = num_x_bins;
	y_step_size = (y_range[1] - y_range[0])/(double)(num_y_bins-1);	
	histogram->y_bin_width = y_step_size; histogram->num_y_bins = num_y_bins;

	//===Fill X Bins===//
	for (i=0; i<histogram->num_x_bins; i++){
		histogram->x_bins[i] = x_range[0] + i*x_step_size;
	}
	histogram->x_bin_range[0] = x_range[0];
	histogram->x_bin_range[1] = x_range[0] + (histogram->num_x_bins)*x_step_size;

	//===Fill Y Bins===//
	for (i=0; i<histogram->num_y_bins; i++){
		histogram->y_bins[i] = y_range[0] + i*y_step_size;
	}
	histogram->y_bin_range[0] = y_range[0];
	histogram->y_bin_range[1] = y_range[0] + (histogram->num_y_bins)*y_step_size;

	//===Initialize Data===//
	initialize_array_uint(histogram->counts, histogram->num_x_bins*histogram->num_y_bins);
	initialize_array_dbl(histogram->probabilities, histogram->num_x_bins*histogram->num_y_bins);
	
	return histogram;
}

void add_to_histogram_2d( histogram_2d_t* self,
						  double complex* points,
						  unsigned int num_points )
{
	unsigned int p, s, real_index, imag_index;
	double distance, min_distance;
	double complex point;

	for (p=0; p<num_points; p++){

		//===Get Data Point===//
		point = points[p];

		//===Check Real Part===//
		distance = 0;
		real_index = UINT_MAX;
		if (creal(point) >= self->x_bin_range[0] && creal(point) <= self->x_bin_range[1]){
			min_distance = DBL_MAX;
			for (s=0; s<self->num_x_bins; s++){
				distance = pow(creal(point) - self->x_bins[s], 2.0);			
				if (distance < min_distance){
					min_distance = distance;
					real_index = s;
				}
			}	
		}

		//===Check Imag Part===//
		distance = 0;
		imag_index = UINT_MAX;
		if (cimag(point) >= self->y_bin_range[0] && cimag(point) <= self->y_bin_range[1]){
			min_distance = DBL_MAX;
			for (s=0; s<self->num_y_bins; s++){
				distance = pow(cimag(point) - self->y_bins[s], 2.0);			
				if (distance < min_distance){
					min_distance = distance;
					imag_index = s;
				}
			}
		}

		//===Incremement Bin===//							 
		if (real_index != UINT_MAX && imag_index != UINT_MAX){
			if (self->counts[real_index + imag_index*self->num_x_bins] < UINT_MAX-1){
				self->counts[real_index + imag_index*self->num_x_bins] += 1;
			}
		}

	}

	return;
}

void make_histogram_binary_2d( histogram_2d_t* self,
							   unsigned int threshold )
{
	
	unsigned int i;
	for (i=0; i<self->num_x_bins*self->num_y_bins; i++){
		if (self->counts[i] >= threshold){
			self->counts[i] = 1;
		}
		else{
			self->counts[i] = 0;
		}
	}

	return;
}

//================================================================================================//
//=================================HISTOGRAM 1D METHODS===========================================//
//================================================================================================//

void initialize_histogram_1d( histogram_1d_t* self,
							  unsigned int num_bins,
							  double* bin_range )
{
	unsigned int i;
	double step_size;

	//===Set Bin Widths===//
	if (num_bins % 2 == 0){
		num_bins++;		
	}
	step_size = (bin_range[1] - bin_range[0])/(double)(num_bins-1);	
	self->bin_width = step_size; self->num_bins = num_bins;

	//===Fill Bins===//
	for (i=0; i<self->num_bins; i++){
		self->bins[i] = bin_range[0] + i*step_size;
	}
	self->start_bin = bin_range[0];
	self->end_bin = bin_range[0] + (self->num_bins)*step_size;
	
	//===Initialize Data===//
	initialize_array_uint(self->counts, self->num_bins);
	initialize_array_dbl(self->probabilities, self->num_bins);
	initialize_array_dbl(self->cdf, self->num_bins);
	
	return;
}

void initialize_histogram_1d_data( histogram_1d_t* self,
						   	  	   unsigned int num_bins,
						   	  	   double* data,
						   	  	   unsigned int length )
{

	unsigned int i, j, closest_bin_index;
	double min, max, distance, min_distance;

	if (data == NULL){
		fprintf(stderr, "Error:: Data Is NULL! In Function -- initialize_histogram_1d_data\n");
		return;
	}

	//===Set Local Data===//
	self->num_bins = num_bins;
	min = find_minimum_dbl(data, length);
	self->start_bin = min;
	max = find_maximum_dbl(data, length);
	self->end_bin = max;

	//===Initialize Probabilities===//
	initialize_array_dbl(self->probabilities, num_bins);

	//===Calculate Bin Width===//
	self->bin_width = (self->end_bin - self->start_bin)/(double)(num_bins);

	//===Create Bins===//	
	for (i=0; i<self->num_bins; i++){
		self->bins[i] = self->start_bin + i*self->bin_width;
	}

	if (data != NULL){
		//===Calculate Bin Counts===//
		for (i=0; i<length; i++){
			//===Find Closest Bin Center===//
			closest_bin_index = 0;
			min_distance = DBL_MAX;
			for(j=0; j<num_bins; j++){
				distance = fabs(data[i] - self->bins[j]); 
				if (distance < min_distance){
					min_distance = distance;
					closest_bin_index = j;
				}
			}
			if (self->counts[closest_bin_index] < UINT_MAX){ 	
				self->counts[closest_bin_index] += 1.0;
			}
			if (self->probabilities[closest_bin_index] < DBL_MAX){ 
				self->probabilities[closest_bin_index] += 1.0;
			}
		}
		//===Normalize Sum To One===//
		normalize_sum_dbl(self->probabilities, self->num_bins);
	}

	return;
}

//Kernel density estimation
//1. for each point xi:
//	 	for many samples (s):
//			Draw a value from the kernel associated with the point. 
//				-- In this case, draw from the Gaussian N(xi,h) centered at xi and of variance h (the bandwidth)
//2. Divide all values in the new array by n -- the sample size
//3. Create a histogram of this new array
//NOTE: Pick the bandwidth to be:
// 		-- h = (4 * sample_std(data)**5 / (3.0 * sample_size) )**1.0/5.0

void initialize_kernel_density_histogram_1d_data( histogram_1d_t* self,
						   				  	 	  unsigned int num_bins,
						   				  	 	  double* data,
						   				  	 	  unsigned int length )
{
	unsigned int i, j, closest_bin_index;
	double min, max, distance, min_distance, std, bandwidth;
	double* kernel_data;
	random_number_generator_t rng;

	//===Error Check===//
	if (data == NULL){
		fprintf(stderr, "Error:: Data Is NULL! In Function -- initialize_kernel_density_histogram\n");
		return;
	}

	//===Make Space For Kernel Data===//
	kernel_data = malloc(length * MAX_KERNEL_SAMPLES * sizeof(*kernel_data));

	//===Make RNG===//
	initialize_random_number_generator(&rng);

	//===Initialize Probabilities===//
	initialize_array_dbl(self->probabilities, num_bins);

	//===Calculate Kernel Bandwidth===//
	std = sqrt(compute_variance_dbl(data,length));
	bandwidth = 4.0 * pow(std, 5.0);
	bandwidth /= (3.0 * (double)length);
	bandwidth = pow(bandwidth, 0.2);

	//===Run Through Data===//
	for (i=0; i<length; i++){
		//===Generate Kernel Samples===//
		for (j=0; j<MAX_KERNEL_SAMPLES; j++){
			kernel_data[j + i*MAX_KERNEL_SAMPLES] = get_random_gaussian(data[i],bandwidth,&rng);
		}
	}

	//===Set Local Data===//
	self->num_bins = num_bins;
	min = find_minimum_dbl(kernel_data, length*MAX_KERNEL_SAMPLES);
	self->start_bin = min;
	max = find_maximum_dbl(kernel_data, length*MAX_KERNEL_SAMPLES);
	self->end_bin = max;

	//===Calculate Bin Width===//
	self->bin_width = (self->end_bin - self->start_bin)/(double)(num_bins);

	//===Create Bins===//	
	for (i=0; i<num_bins; i++){
		self->bins[i] = self->start_bin + i*self->bin_width;
	}

	//===Run Through Kernel Data===//
	for (i=0; i<length*MAX_KERNEL_SAMPLES; i++){
		//===Find Closest Bin Center===//
		closest_bin_index = 0;
		min_distance = DBL_MAX;
		for(j=0; j<num_bins; j++){
			distance = fabs(kernel_data[i] - self->bins[j]); 
			if (distance < min_distance){
				min_distance = distance;
				closest_bin_index = j;
			}
		}
		if (self->counts[closest_bin_index] < UINT_MAX){
			self->counts[closest_bin_index] += 1.0;
		}
		if (self->probabilities[closest_bin_index] < DBL_MAX){
			self->probabilities[closest_bin_index] += 1.0;
		}
	}

	//===Normalize Sum To One===//
	normalize_sum_dbl(self->probabilities, self->num_bins);

	//===Clean Up===//
	free(kernel_data);

	return;
}

void add_to_histogram_1d( histogram_1d_t* self,
						  double* data,
						  unsigned int data_length )
{
	unsigned int i, j, closest_bin_index;
	double min_distance, distance;

	if (data != NULL){
		//===Calculate Bin Counts===//
		for (i=0; i<data_length; i++){
			//===Find Closest Bin Center===//
			closest_bin_index = 0;
			min_distance = DBL_MAX;
			for(j=0; j<self->num_bins; j++){
				distance = fabs(data[i] - self->bins[j]); 
				if (distance < min_distance){
					min_distance = distance;
					closest_bin_index = j;
				}
			}
			if (self->counts[closest_bin_index] < UINT_MAX){ 	
				self->counts[closest_bin_index] += 1.0;
			}
		}
	}

	return;
}

//this one adds the data to the counts, but does so with a kde
void add_to_histogram_1d_kde( histogram_1d_t* self,
							  double* data,
							  unsigned int data_length )
{
	unsigned int i, j, closest_bin_index;
	double std, bandwidth, min_distance, distance;
	double* kernel_data;
	random_number_generator_t rng;

	//===Make Space For Kernel Data===//
	kernel_data = malloc(data_length * MAX_KERNEL_SAMPLES * sizeof(*kernel_data));

	//===Make RNG===//
	initialize_random_number_generator(&rng);

	//===Calculate Kernel Bandwidth===//
	std = sqrt(compute_variance_dbl(data,data_length));
	bandwidth = 4.0 * pow(std, 5.0);
	bandwidth /= (3.0 * (double)data_length);
	bandwidth = pow(bandwidth, 0.2);

	//===Run Through Data===//
	for (i=0; i<data_length; i++){
		//===Generate Kernel Samples===//
		for (j=0; j<MAX_KERNEL_SAMPLES; j++){
			kernel_data[j + i*MAX_KERNEL_SAMPLES] = get_random_gaussian(data[i],bandwidth,&rng);
		}
	}

	//===Run Through Kernel Data===//
	for (i=0; i<data_length*MAX_KERNEL_SAMPLES; i++){
		//===Find Closest Bin Center===//
		closest_bin_index = 0;
		min_distance = DBL_MAX;
		for(j=0; j<self->num_bins; j++){
			distance = fabs(kernel_data[i] - self->bins[j]); 
			if (distance < min_distance){
				min_distance = distance;
				closest_bin_index = j;
			}
		}
		if (self->counts[closest_bin_index] < UINT_MAX){
			self->counts[closest_bin_index] += 1.0;
		}
	}

	//===Clean Up===//
	free(kernel_data);

	return;
}

void initialize_histogram_probabilities( histogram_1d_t* self )
{
	unsigned int i;
	for (i=0; i<self->num_bins; i++) self->probabilities[i] = self->counts[i];
	normalize_sum_dbl(self->probabilities, self->num_bins);
	return;
}

void print_histogram_1d( histogram_1d_t* self,
						 char* filename )
{
	unsigned int i;
	FILE *fout;

	//===Sanity Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Histogram Is NULL! In Function -- print_histogram_1d!\n");
		quit();
	}
	//===Print Shit===//
	fout = fopen(filename,"w");
	for (i=0; i<self->num_bins; i++){
		fprintf(fout, "%lf %lf %d\n", self->bins[i], self->probabilities[i], self->counts[i]);
	}
	fclose(fout);
	return;
}

double compute_histogram_1d_entropy( histogram_1d_t* histogram )
{
	unsigned int i;
	double temp, entropy;
	entropy = 0;	
	for (i=0; i<histogram->num_bins; i++){
		if (fabs(histogram->probabilities[i]) > EPS){
			temp = (-1.0 * histogram->probabilities[i] * log2(histogram->probabilities[i]));
			entropy += temp; 
		}
	}
	return entropy;
}

double compute_histogram_1d_cross_entropy( histogram_1d_t* histogram1, 
										   histogram_1d_t* histogram2 )
{
	unsigned int i, j, closest_bin_index;
	double distance, min_distance, cross_entropy;

	cross_entropy = 0;
	for (i=0; i<MIN(histogram1->num_bins, histogram2->num_bins); i++){

		//===Find Closest Bin Of Histogram 2===//
		closest_bin_index = 0;
		min_distance = DBL_MAX;
		for (j=0; j<histogram2->num_bins; j++){
			distance = fabs( histogram1->bins[i] - histogram2->bins[j] );
			if (distance < min_distance){
				min_distance = distance;
				closest_bin_index = j;
			}
		}

		//===Compute Cross Entropy===//
		if (fabs(histogram1->probabilities[i]) > 0 && fabs(histogram2->probabilities[closest_bin_index] > 0)){
			cross_entropy += -1.0 * histogram1->probabilities[i] * log2(histogram2->probabilities[closest_bin_index]);
		}
	}

	return cross_entropy;
}


//Sample Kullback-Leibler Divergence D(P || Q) -- divergence of approximated Q from true P
// Q -- approximated arma model && P -- true data set
// Create the histogram of the data sets P & Q -- maybe even create a KDE of the data as the histogram
// then calculate:
//		1. H(P): Entropy of P -- sum_{x} p(x) * log(p(x)) -- where x is over bins of histogram, and p(x) is value of histogram
//		2. H(P,Q): Cross Entropy of P and Q -- sum_{x} p(x) * log(q(x)) -- so make sure histograms have same bin starts/widths/ends!!! 

double compute_sample_kullback_leibler_divergence_1d( histogram_1d_t* true_histogram,
												   	  histogram_1d_t* approx_histogram )
{
	double divergence;
	
	divergence = compute_histogram_1d_cross_entropy(true_histogram, approx_histogram);
	divergence -= compute_histogram_1d_entropy(true_histogram);
	
	return divergence;
}

double compute_jensen_shannon_divergence_1d( histogram_1d_t* true_histogram,
										     histogram_1d_t* approx_histogram )
{
	unsigned int i;
	double divergence;
	histogram_1d_t temp;

	//===Compute Smoothed Histogram===//
	for (i=0; i<approx_histogram->num_bins; i++){
		temp.counts[i] = (true_histogram->counts[i] + approx_histogram->counts[i])/2;
		temp.cdf[i] = (true_histogram->cdf[i] + approx_histogram->cdf[i])/2;
		temp.probabilities[i] = (true_histogram->probabilities[i] + approx_histogram->probabilities[i])/2;
	}

	//===Compute JS Divergence===//
	divergence = compute_sample_kullback_leibler_divergence_1d(true_histogram, &temp);
	divergence += compute_sample_kullback_leibler_divergence_1d(approx_histogram, &temp);
	return divergence;
}

void compute_histogram_cdf( histogram_1d_t* self )
{
	unsigned int i;
	double sum;

	//===Initialize Probabilities===//
	initialize_histogram_probabilities(self);

	//===Sum Up Probabilities===//
	sum = 0;
	for (i=0; i<self->num_bins; i++){ 
		sum += self->probabilities[i];
		self->cdf[i] = sum;
	}
		

	return;
}

//Correlation Coefficient is:
// cos(theta) = <x,y>/(<x,x> * <y,y>)

double compute_correlation_coefficient_dbl( double* signal1,
											double* signal2,
											unsigned int signal_length )
{
	double xy, xx, yy;

	//===Compute Covariance===//
	xy = compute_inner_product_dbl(signal1, signal2, signal_length);

	//===Compute Variances===//
	xx = sqrt(compute_inner_product_dbl(signal1, signal1, signal_length));
	yy = sqrt(compute_inner_product_dbl(signal2, signal2, signal_length));
	
	return xy/(xx * yy);
}

double compute_jaccard_index_dbl( double* signal1,
								  double* signal2,
								  unsigned int signal_length )
{
	double xy, xx, yy;

	//===Compute Covariance===//
	xy = compute_inner_product_dbl(signal1, signal2, signal_length);

	//===Compute Variances===//
	xx = (compute_inner_product_dbl(signal1, signal1, signal_length));
	yy = (compute_inner_product_dbl(signal2, signal2, signal_length));

	return xy/(xx + yy - xy);
}

double compute_sorenson_index_dbl( double* signal1,
								   double* signal2,
								   unsigned int signal_length )
{
	double xy, xx, yy;

	//===Compute Covariance===//
	xy = compute_inner_product_dbl(signal1, signal2, signal_length);

	//===Compute Variances===//
	xx = (compute_inner_product_dbl(signal1, signal1, signal_length));
	yy = (compute_inner_product_dbl(signal2, signal2, signal_length));

	return (2.0*xy)/(xx + yy);
}


double compute_anderson_darling_statistic_dbl( histogram_1d_t* empirical_histogram,
										   	   histogram_1d_t* true_histogram )
{
	unsigned int i;
	double temp, stat;
	stat = 0;
	for (i=1; i<empirical_histogram->num_bins; i++){
		temp = pow( (empirical_histogram->cdf[i] - true_histogram->cdf[i]) + EPS, 2.0) + EPS;
		temp /= (((true_histogram->cdf[i] +EPS) * (1.0 - (true_histogram->cdf[i] + EPS))) + EPS);
		stat += temp;
	}
	stat *= (double)(empirical_histogram->num_bins);
	return stat;
}

double compute_kuiper_statistic_dbl( histogram_1d_t* empirical_histogram,
									 histogram_1d_t* true_histogram )
{
	double dplus, dminus;
	double temp_result[MAX_NUM_BINS];

	//===Find Maximum Where Empirical > True===//
	subtract_vectors_dbl(empirical_histogram->cdf, true_histogram->cdf, 
						 empirical_histogram->num_bins, temp_result);
	dplus = find_maximum_dbl(temp_result, empirical_histogram->num_bins);

	//===Find Maximum Where True > Empirical===//
	subtract_vectors_dbl(true_histogram->cdf, empirical_histogram->cdf, 
						 true_histogram->num_bins, temp_result);
	dminus = find_maximum_dbl(temp_result, true_histogram->num_bins);

	return dplus + dminus;
}


double compute_hellinger_statistic_dbl(	histogram_1d_t* empirical_histogram,
								 	 	histogram_1d_t* true_histogram )
{
	unsigned int i;
	double hellinger_distance;
	hellinger_distance = 0;
	for (i=0; i<empirical_histogram->num_bins; i++){
		hellinger_distance += (sqrt(empirical_histogram->probabilities[i]) - sqrt(true_histogram->probabilities[i]));
	}
	hellinger_distance = sqrt(hellinger_distance);
	hellinger_distance /= sqrt(2.0);
	return hellinger_distance;
}

double compute_total_variation_statistic_dbl( histogram_1d_t* empirical_histogram,
								 	 		  histogram_1d_t* true_histogram )
{
	
	unsigned int i;
	double total_variation_distance;
	total_variation_distance = 0;
	for (i=0; i<empirical_histogram->num_bins; i++){
		total_variation_distance += fabs(empirical_histogram->probabilities[i] - true_histogram->probabilities[i]);
	}
	total_variation_distance /= 2.0;
	return total_variation_distance;
}

double compute_separation_statistic_dbl( histogram_1d_t* empirical_histogram,
								 	 	 histogram_1d_t* true_histogram )
{
	
	unsigned int i;
	double separation, max_separation;
	max_separation = -(DBL_MAX-1000);
	for (i=0; i<empirical_histogram->num_bins; i++){
		separation = 1.0 - (empirical_histogram->probabilities[i]+EPS)/(true_histogram->probabilities[i]+EPS);
		if (separation > max_separation){
			max_separation = separation;
		}
	}
	return max_separation;
}

double compute_chi_squared_statistic_dbl( histogram_1d_t* empirical_histogram,
										  histogram_1d_t* true_histogram )
{

	unsigned int i;
	double temp, chi_squared_distance;
	chi_squared_distance = 0;
	for (i=0; i<empirical_histogram->num_bins; i++){
		temp = pow(empirical_histogram->probabilities[i] - true_histogram->probabilities[i], 2.0) + EPS;
		temp /= (true_histogram->probabilities[i] + EPS);
		chi_squared_distance += temp;
	}
	return chi_squared_distance;
}

double compute_bhattacharyya_distance_dbl( histogram_1d_t* empirical_histogram,
										   histogram_1d_t* true_histogram )
{

	unsigned int i;
	double temp, bhattacharyya_distance;
	bhattacharyya_distance = 0;
	for (i=0; i<empirical_histogram->num_bins; i++){
		temp = sqrt(empirical_histogram->probabilities[i] * true_histogram->probabilities[i]) + EPS;
		bhattacharyya_distance += temp;
	}
	return bhattacharyya_distance;
}

double compute_bhattacharyya_angle_dbl( histogram_1d_t* empirical_histogram,
										histogram_1d_t* true_histogram )
{

	unsigned int i;
	double temp, bhattacharyya_distance;
	bhattacharyya_distance = 0;
	for (i=0; i<empirical_histogram->num_bins; i++){
		temp = sqrt(empirical_histogram->probabilities[i] * true_histogram->probabilities[i]) + EPS;
		bhattacharyya_distance += temp;
	}
	return acos(bhattacharyya_distance);
}

double compute_kantorovich_metric_dbl(histogram_1d_t* empirical_histogram,
									  histogram_1d_t* true_histogram )
{

	unsigned int i;
	double kantorovich_metric;
	kantorovich_metric = 0;
	for (i=0; i<empirical_histogram->num_bins; i++){
		kantorovich_metric += fabs(empirical_histogram->cdf[i] - true_histogram->cdf[i]);
	}
	return kantorovich_metric;
}
