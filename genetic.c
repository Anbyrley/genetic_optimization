#include "genetic.h"



//================================================================================================//
//=========================================GENE METHODS===========================================//
//================================================================================================//

void initialize_gene( gene_t* self,
					  double* bounds,
					  double allele,
					  double mutation_rate,
					  double mutation_magnitude )
{
	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Gene Is NULL! In Function -- initialize_gene!\n");
		quit();
	}
	
	//===Set Locals===//
	self->bounds[0] = bounds[0]; self->bounds[1] = bounds[1];
	self->allele = allele;
	self->mutation_rate = mutation_rate;
	self->mutation_magnitude = mutation_magnitude;

	return;
}

void copy_gene( gene_t* original,
				gene_t* copy )
{
	copy->mutation_rate = original->mutation_rate;
	copy->mutation_magnitude = original->mutation_magnitude;
	copy_array_dbl(original->bounds, 2, copy->bounds);
	copy->allele = original->allele;

	return;
}

void mutate_gene( gene_t* self,
				  random_number_generator_t* rng )
{

	double random_number;
	double domain_center;

	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Gene Is NULL! In Function -- mutate_gene!\n");
		quit();
	}
	if (rng == NULL){
		fprintf(stderr, "Error:: Random Number Generator Is NULL! In Function -- mutate_gene!\n");
		quit();
	}

	//===Generate Mutation===//
	#if ALLOW_GAUSSIAN_MUTATION
		#if ALLOW_DOMAIN_CENTERED_MUTATION
			domain_center = 0.5*(self->bounds[0] + self->bounds[1]);
		#else
			domain_center = self->allele;
		#endif
		random_number = get_random_gaussian(domain_center, self->mutation_magnitude, rng);
		//===Set Mutated Feature===//
		self->allele = MAX((MIN(self->bounds[1], random_number)),self->bounds[0]);
	#else
		#if ALLOW_DOMAIN_CENTERED_MUTATION
			domain_center = self->bounds[0];
			random_number = get_random_uniform_dbl(domain_center, self->bounds[1], rng);
			//===Set Mutated Feature===//
			self->allele = MAX((MIN(self->bounds[1], random_number)),self->bounds[0]);			
		#else
			domain_center = MIN((self->allele-self->bounds[0]), (self->bounds[1]-self->allele));
			random_number = get_random_uniform_dbl(self->allele-fabs(domain_center), self->allele+fabs(domain_center), rng);
			//===Set Mutated Feature===//
			self->allele = MAX((MIN(self->bounds[1], random_number)),self->bounds[0]);
		#endif
	#endif

	return;
}



//================================================================================================//
//======================================CHROMOSOME METHODS========================================//
//================================================================================================//

void initialize_chromosome( chromosome_t* self,
							double* gene_bounds,
							double* mutation_magnitudes,
							double* mutation_rates,
							unsigned int num_genes,
							random_number_generator_t* rng )
{
	unsigned int g;
	double allele, mutation_magnitude, mutation_rate;

	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Chromosome Is NULL! In Function -- initialize_chromosome!\n");
		quit();
	}

	//===Set Locals===//
	self->num_genes = num_genes;
	self->fitness = 0;

	//===Set Dominance===//
	self->dominant = 2;
	#if ALLOW_DOMINANT_GENETICS
		self->dominant = get_random_uniform(0, 100, rng);
		if (self->dominant >= 50){
			self->dominant = 1;
		}
		else{
			self->dominant = 0;
		}
	#endif
	
	//===Initialize Genes===//
	for (g=0; g<self->num_genes; g++){

		//===Set Mutation Information===//
		if (mutation_magnitudes == NULL){
		#if ALLOW_MUTATION_MAGNITUDE_RECOMBINATION
			mutation_magnitude = get_random_uniform_dbl(0,1,rng)*(fabs(self->genes[g].bounds[1] - self->genes[g].bounds[0]));
		#else
			mutation_magnitude = 0.1*(fabs(self->genes[g].bounds[1] - self->genes[g].bounds[0]));
		#endif
		}
		else{
			mutation_magnitude = mutation_magnitudes[g];
		}
		if (mutation_rates == NULL){
		#if ALLOW_MUTATION_RATE_RECOMBINATION
			mutation_rate = get_random_uniform_dbl(0,1,rng);
		#else
			mutation_rate = 0.2;
		#endif
		}
		else{
			mutation_rate = mutation_rates[g];
			if (mutation_rate > 1){
				fprintf(stderr, "Error:: Mutation Rate Is Greater Than 1! In Function -- initialize_chromosome!\n");
				quit();
			}
		}

		//===Generate Gene's Allele===//
		//need to think about how to get gene bounds into the problem//
		allele = get_random_uniform_dbl(gene_bounds[2*g+0], gene_bounds[2*g+1], rng);

		//===Initialize Chromosome's Gene===//
		initialize_gene(&(self->genes[g]), (gene_bounds + 2*g), allele, mutation_rate, mutation_magnitude);
	}

	return;
}

void print_chromosome( chromosome_t* self,
					   FILE* file )
{
	unsigned int g;
	if (isatty(fileno(file))){
		//===Print To Stdout===//
		fprintf(file, "Genes:	");	
		for (g=0; g<self->num_genes; g++){
			fprintf(file, "%+10.5f ", self->genes[g].allele);
		}
		fprintf(file, "\nFitness: %+10.5f ", self->fitness);	
		fprintf(file, "\n");
	}
	else{
		//===Print To File===//
		fprintf(file, "Genes:	");	
		for (g=0; g<self->num_genes; g++){
			fprintf(file, "%+10.5f ", self->genes[g].allele);
		}
		fprintf(file, "  Fitness: %+10.5f ", self->fitness);	
		fprintf(file, "\n");
	}
	return;	
}

void copy_chromosome( chromosome_t* original,
					  chromosome_t* copy )
{
	unsigned int g;

	copy->num_genes = original->num_genes;
	copy->dominant = original->dominant;
	copy->fitness = original->fitness;
	
	for (g=0; g<original->num_genes; g++){
		copy_gene(&(original->genes[g]), &(copy->genes[g]));
	}

	#if TEST_GENETIC_OPTIMIZER
		copy->calculate_fitness = original->calculate_fitness;
	#endif

	return;
}


void mutate_chromosome_genes( chromosome_t* self,
							  random_number_generator_t* rng )
{

	unsigned int g;
	double random_number;

	for (g=0; g<self->num_genes; g++){	

		//===Generate Random Number===//
		random_number = get_random_uniform_dbl(0, 1, rng);

		//===Check If Need To Mutate===//
		if (random_number <= self->genes[g].mutation_rate){
			mutate_gene(&(self->genes[g]), rng);
		}
	}

	//===Recalculate Fitness===//
	#if TEST_GENETIC_OPTIMIZER
		self->calculate_fitness(self);
	#endif

	return;
}

unsigned int are_chromosomes_species_compatibile( chromosome_t* chromosome_1,
												  chromosome_t* chromosome_2 )
{
	if (chromosome_1 == NULL){
		fprintf(stderr, "Error:: Chromosome_1 Is NULL 1! In Function -- are_chromosomes_species_compatibile!\n");
		quit();
	}
	if (chromosome_2 == NULL){
		fprintf(stderr, "Error:: Chromosome_2 Is NULL 1! In Function -- are_chromosomes_species_compatibile!\n");
		quit();
	}
	return 1;
}

unsigned int are_chromosomes_niches_compatibile( chromosome_t* chromosome_1,
												 chromosome_t* chromosome_2 )
{
	if (chromosome_1 == NULL){
		fprintf(stderr, "Error:: Chromosome_1 Is NULL 1! In Function -- are_chromosomes_species_compatibile!\n");
		quit();
	}
	if (chromosome_2 == NULL){
		fprintf(stderr, "Error:: Chromosome_2 Is NULL 1! In Function -- are_chromosomes_species_compatibile!\n");
		quit();
	}
	return 1;
}

unsigned int chromosomes_are_equal( chromosome_t* chromosome1,
									chromosome_t* chromosome2 )
{
	unsigned int g;

	//===Error Check===//
	if (chromosome1 == NULL){
		fprintf(stderr, "Error:: Chromosome1 Is NULL! In Function -- chromosomes_are_equal!\n");
		quit();
	}
	if (chromosome2 == NULL){
		fprintf(stderr, "Error:: Chromosome2 Is NULL! In Function -- chromosomes_are_equal!\n");
		quit();
	}

	//===Run Comparison===//
	if (chromosome1->num_genes != chromosome2->num_genes){
		return 0;
	}
	for (g=0; g<chromosome1->num_genes; g++){
		if (fabs(chromosome1->genes[g].allele - chromosome2->genes[g].allele) > 2.0*EPS){
			return 0;
		}
	}
	return 1;
}


void breed_chromosomes_fuzzy( chromosome_t* parents,
							  chromosome_t* children,
							  unsigned int num_parents,
							  random_number_generator_t* rng )
{

	unsigned int g;		
	double alpha;
	chromosome_t parent_1, parent_2;
	#if ALLOW_MULTIPLE_PARENTS
		unsigned int parent_indices[2];
	#endif

	//===Error Checks===//
	if (num_parents < 2){
		fprintf(stderr, "Error:: Need At Least Two Parents! In Function -- breed_chromosomes_fuzzy!\n");
		quit();
	}	

	//===Generate Alpha===//
	alpha = get_random_uniform_dbl(-1, 1, rng);
	if (alpha >= -1 && alpha < 0){
		alpha = 1 + alpha;
	}
	else{
		alpha = 1 - alpha;
	}

	//===Copy Over===//
	copy_chromosome(&(parents[0]), &(parent_1));
	copy_chromosome(&(parents[1]), &(parent_2));

	//===Copy Base Info===//
	copy_chromosome(&(parent_1), &(children[0]));
	copy_chromosome(&(parent_2), &(children[1]));
	
	#if !ALLOW_MULTIPLE_PARENTS
		num_parents = 2;
		if (parent_1.num_genes != parent_2.num_genes){
			fprintf(stderr, "Error:: Parents Have An Unequal Number Of Genes! In Function -- breed_chromosomes_fuzzy!\n");
			quit();
		}
		#if !ALLOW_MUTATION_RATE_RECOMBINATION
			for (g=0; g<parent_1.num_genes; g++){
				if (fabs(parent_1.genes[g].mutation_rate - parent_2.genes[g].mutation_rate)>2*EPS 
					&& !ALLOW_MUTATION_RATE_RECOMBINATION){
					fprintf(stderr, "Error:: Parents Genes Have Different Mutation Rates! In Function -- breed_chromosomes_fuzzy!\n");
					quit();
				}
			}
		#endif
		#if !ALLOW_MUTATION_MAGNITUDE_RECOMBINATION
			for (g=0; g<parent_1.num_genes; g++){
				if (fabs(parent_1.genes[g].mutation_magnitude - parent_2.genes[g].mutation_magnitude)>2*EPS 
					&& !ALLOW_MUTATION_MAGNITUDE_RECOMBINATION){
					fprintf(stderr, "Error:: Parents Genes Have Different Mutation Magnitude! In Function -- breed_chromosomes_fuzzy!\n");
					quit();
				}
			}
		#endif
	#endif


	//===Fuzzy Recombo===//
	for (g=0; g<(parent_1.num_genes); g++){

		#if ALLOW_MULTIPLE_PARENTS

			//===Pick Parents===//
			parent_indices[0] = 0;
			parent_indices[1] = 0;
			while(parent_indices[0] == parent_indices[1]){
				parent_indices[0] = get_random_uniform(0, 1000, rng) % num_parents;
				parent_indices[1] = get_random_uniform(0, 1000, rng) % num_parents;
			}
			copy_chromosome(&(parents[parent_indices[0]]), &(parent_1));
			copy_chromosome(&(parents[parent_indices[1]]), &(parent_2));

			//===Error Checks===//
			if (parent_1.num_genes != parent_2.num_genes){
				fprintf(stderr, "Error:: Parents Have An Unequal Number Of Genes! In Function -- breed_chromosomes_fuzzy!\n");
				quit();
			}
			#if !ALLOW_MUTATION_RATE_RECOMBINATION
				for (g=0; g<parent_1.num_genes; g++){
					if (fabs(parent_1.genes[g].mutation_rate - parent_2.genes[g].mutation_rate)>2*EPS 
						&& !ALLOW_MUTATION_RATE_RECOMBINATION){
					   fprintf(stderr, "Error:: Parents Genes Have Different Mutation Rates! In Function -- breed_chromosomes_fuzzy!\n");
						quit();
					}
				}
			#endif
			#if !ALLOW_MUTATION_MAGNITUDE_RECOMBINATION
				for (g=0; g<parent_1.num_genes; g++){
					if (fabs(parent_1.genes[g].mutation_magnitude - parent_2.genes[g].mutation_magnitude)>2*EPS 
						&& !ALLOW_MUTATION_MAGNITUDE_RECOMBINATION){
					  fprintf(stderr, "Error:: Parents Genes Have Different Mutation Magnitude! In Function -- breed_chromosomes_fuzzy!\n");
						quit();
					}
				}
			#endif

		#endif

		//===Recombo Child 1===//
		children[0].genes[g].allele = alpha * parent_1.genes[g].allele;
		children[0].genes[g].allele += (1.0 - alpha) * parent_2.genes[g].allele;

		//===Recombo Child 2===//
		children[1].genes[g].allele = (1.0 - alpha) * parent_1.genes[g].allele;
		children[1].genes[g].allele += alpha * parent_2.genes[g].allele;

	}


	//===Recombo Mutation Rates===//
	#if ALLOW_MUTATION_RATE_RECOMBINATION

		#if ALLOW_MULTIPLE_PARENTS
			//===Pick Parents===//
			parent_indices[0] = 0;
			parent_indices[1] = 0;
			while(parent_indices[0] == parent_indices[1]){
				parent_indices[0] = get_random_uniform(0, 1000, rng) % num_parents;
				parent_indices[1] = get_random_uniform(0, 1000, rng) % num_parents;
			}
			copy_chromosome(&(parents[parent_indices[0]]), &(parent_1));
			copy_chromosome(&(parents[parent_indices[1]]), &(parent_2));
		#endif

		//===Generate Alpha===//
		alpha = get_random_uniform_dbl(-1, 1, rng);
		if (alpha >= -1 && alpha < 0){
			alpha = 1 + alpha;
		}
		else{
			alpha = 1 - alpha;
		}
		for (g=0; g<parent_1.num_genes; g++){

			//===Recombo Child 1===//
			children[0].genes[g].mutation_rate = alpha * parent_1.genes[g].mutation_rate;
			children[0].genes[g].mutation_rate += (1.0-alpha) * parent_2.genes[g].mutation_rate;

			//===Recombo Child 2===//
			children[1].genes[g].mutation_rate = (1.0-alpha) * parent_1.genes[g].mutation_rate;
			children[1].genes[g].mutation_rate += alpha * parent_2.genes[g].mutation_rate;

		}
	#endif

	//===Recombo Mutation Magnitudes===//
	#if ALLOW_MUTATION_MAGNITUDE_RECOMBINATION

		#if ALLOW_MULTIPLE_PARENTS
			//===Pick Parents===//
			parent_indices[0] = 0;
			parent_indices[1] = 0;
			while(parent_indices[0] == parent_indices[1]){
				parent_indices[0] = get_random_uniform(0, 1000, rng) % num_parents;
				parent_indices[1] = get_random_uniform(0, 1000, rng) % num_parents;
			}
			copy_chromosome(&(parents[parent_indices[0]]), &(parent_1));
			copy_chromosome(&(parents[parent_indices[1]]), &(parent_2));
		#endif

		//===Generate Alpha===//
		alpha = get_random_uniform_dbl(-1, 1, rng);
		if (alpha >= -1 && alpha < 0){
			alpha = 1 + alpha;
		}
		else{
			alpha = 1 - alpha;
		}
		for (g=0; g<parent_1.num_genes; g++){

			//===Recombo Child 1===//
			children[0].genes[g].mutation_magnitude = alpha * parent_1.genes[g].mutation_magnitude;
			children[0].genes[g].mutation_magnitude += (1.0-alpha) * parent_2.genes[g].mutation_magnitude;

			//===Recombo Child 2===//
			children[1].genes[g].mutation_magnitude = (1.0-alpha) * parent_1.genes[g].mutation_magnitude;
			children[1].genes[g].mutation_magnitude += alpha * parent_2.genes[g].mutation_magnitude;

		}
	#endif

	#if ALLOW_DOMINANT_GENETICS
		children[0].dominant = get_random_uniform(0, 100, rng);
		if (children[0].dominant >= 50){
			children[0].dominant = 1;
		}
		else{
			children[0].dominant = 0;
		}

		children[1].dominant = get_random_uniform(0, 100, rng);
		if (children[1].dominant >= 50){
			children[1].dominant = 1;
		}
		else{
			children[1].dominant = 0;
		}
	#endif

	//===Calculate Fitness===//
	#if TEST_GENETIC_OPTIMIZER
		children[0].calculate_fitness(&(children[0]));
		children[1].calculate_fitness(&(children[1]));
	#else
		children[0].fitness = 0;
		children[1].fitness = 0;
	#endif
	
	return;
}


void sort_chromosomes_by_fitness( chromosome_t* chromosomes,
							   unsigned int num_chromosomes )
{
	unsigned int c;
	unsigned int sorted_fitness_indices[(MAX_NUM_CHROMOSOMES+1)];
	double chromosomes_fitness[(MAX_NUM_CHROMOSOMES+1)];
	chromosome_t temp_chromosomes[(MAX_NUM_CHROMOSOMES+1)];

	//===Copy All Chromosomes===//
	for (c=0; c<num_chromosomes; c++){
		copy_chromosome(&(chromosomes[c]), &(temp_chromosomes[c]));
	}
	
	//===Get All Fitness Values===//
	for (c=0; c<num_chromosomes; c++){
		chromosomes_fitness[c] = chromosomes[c].fitness;
	}

	//===Sort Fitness Values In Decreasing Order===//
	sort_array_indices_dbl(chromosomes_fitness, num_chromosomes, sorted_fitness_indices);
	reverse_array_uint(sorted_fitness_indices, num_chromosomes);

	//===Copy Back Chromosomes In Sorted Order===//
	for (c=0; c<num_chromosomes; c++){
		copy_chromosome(&(temp_chromosomes[sorted_fitness_indices[c]]), &(chromosomes[c]));
	}

	return;
}


void shuffle_chromosomes( chromosome_t* chromosomes,
						  unsigned int num_chromosomes,
						  random_number_generator_t* rng )
{
	unsigned int c;
	unsigned int random_indices[(MAX_NUM_CHROMOSOMES+1)];
	chromosome_t temp_chromosomes[(MAX_NUM_CHROMOSOMES+1)];

	//===Copy Chromosomes===//
	for (c=0; c<num_chromosomes; c++){
		copy_chromosome(&(chromosomes[c]), &(temp_chromosomes[c]));
	}

	//===Get Indices===//
	for (c=0; c<num_chromosomes; c++) random_indices[c] = c;
	shuffle_array_uint(random_indices, num_chromosomes, rng);

	//===Copy Back===//
	for (c=0; c<num_chromosomes; c++){
		copy_chromosome(&(temp_chromosomes[random_indices[c]]), &(chromosomes[c]));
	}
	
	return;
}

void find_most_fit_chromosome( chromosome_t* chromosomes,
							   unsigned int num_chromosomes,
							   chromosome_t* most_fit_chromosome,
							   unsigned int* most_fit_chromosome_index )
{
	unsigned int i, max_fitness_index;
	double max_fitness;

	//===Find Max Fitness===//
	max_fitness = -(DBL_MAX-1);
	for (i=0; i<num_chromosomes; i++){
		if (chromosomes[i].fitness > max_fitness){
			max_fitness = chromosomes[i].fitness;
			max_fitness_index = i;
		}
	}

	//===Copy Over===//
	if (most_fit_chromosome != NULL){
		copy_chromosome(&(chromosomes[max_fitness_index]), most_fit_chromosome);
	}	
	if (most_fit_chromosome_index != NULL){
		*most_fit_chromosome_index = max_fitness_index;
	}

	return;
}

void find_least_fit_chromosome( chromosome_t* chromosomes,
							    unsigned int num_chromosomes,
							    chromosome_t* least_fit_chromosome,
							    unsigned int* least_fit_chromosome_index )
{
	unsigned int i, min_fitness_index;
	double min_fitness;

	//===Find Min Fitness===//
	min_fitness = (DBL_MAX-1);
	for (i=0; i<num_chromosomes; i++){
		if (chromosomes[i].fitness < min_fitness){
			min_fitness = chromosomes[i].fitness;
			min_fitness_index = i;
		}
	}

	//===Copy Over===//
	if (least_fit_chromosome != NULL){
		copy_chromosome(&(chromosomes[min_fitness_index]), least_fit_chromosome);
	}	
	if (least_fit_chromosome_index != NULL){
		*least_fit_chromosome_index = min_fitness_index;
	}

	return;
}


//================================================================================================//
//======================================OPTIMIZER METHODS=========================================//
//================================================================================================//

void initialize_genetic_optimizer( genetic_optimizer_t* self,
								   double* gene_bounds, 
								   double* mutation_magnitudes, 
								   double* mutation_rates,
								   unsigned int num_genes,
								   unsigned int num_chromosomes,
								   unsigned int num_elites,
								   unsigned int num_parents,
								   unsigned int generation_size )
{
	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Genetic Optimizer Is NULL! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}
	if (gene_bounds == NULL){
		fprintf(stderr, "Error:: Gene Bounds Are NULL! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}
	if (num_genes > MAX_NUM_GENES){
		fprintf(stderr, "Error:: Too Many Genes! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}
	if (num_chromosomes > (MAX_NUM_CHROMOSOMES+1)){
		fprintf(stderr, "Error:: Too Many Chromosomes! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}
	if (num_elites >= num_chromosomes){
		fprintf(stderr, "Error:: Too Many Elites! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}
	if (generation_size > num_chromosomes || generation_size == 0){
		fprintf(stderr, "Error:: Bad Generation Size! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}
	if (num_parents < 2){
		fprintf(stderr, "Error:: Must Have At Least Two Parents! In Function -- initialize_genetic_optimizer!\n");
		quit();
	}

	unsigned int c;

	//===Set Locals===//
	self->num_chromosomes = num_chromosomes;
	self->num_elites = num_elites;
	self->generation_size = generation_size;
	self->num_parents = num_parents;

	//===Initialize RNG===//
	initialize_random_number_generator(&(self->rng));

	//===Initialize Chromosomes===//
	for (c=0; c<num_chromosomes; c++){ 
		initialize_chromosome(&(self->chromosomes[c]), gene_bounds, mutation_magnitudes, mutation_rates, num_genes,  &(self->rng));
	}
	initialize_chromosome(&(self->global_best_chromosome), gene_bounds, mutation_magnitudes, mutation_rates, num_genes,  &(self->rng));
	self->global_best_chromosome.fitness = -(DBL_MAX-1);		

	//===Initialize Fitness History===//
	initialize_ring_buffer(&(self->fitness_history), 50);

	return;
}

void print_genetic_population( genetic_optimizer_t* self,
							   FILE* file )
{
	unsigned int c;	
	if (isatty(fileno(file))){
		//===Print To Stdout===//	
		for (c=0; c<self->num_chromosomes; c++){
			fprintf(file, "Chromosome %d:	",c);
			print_chromosome(&(self->chromosomes[c]), file);
			//fprintf(file, "\n");
		}
	}
	else{
		//===Print To File===//
		for (c=0; c<self->num_chromosomes; c++){
			print_chromosome(&(self->chromosomes[c]), file);
		}
	}


	return;
}


void scale_genetic_population_fitness( genetic_optimizer_t* self )
{
	unsigned int c;
	double fitness[(MAX_NUM_CHROMOSOMES+1)];
	double temp_dbl, fitness_mean, fitness_std;

	if (self == NULL){
		fprintf(stderr, "Error:: Genetic Optimizer Is NULL! In Function -- scale_genetic_population_fitness!\n");
		quit();
	}

	//===Get All Fitness Values===//
	for (c=0; c<self->num_chromosomes; c++){
		fitness[c] = self->chromosomes[c].fitness;
	}
	fitness[self->num_chromosomes] = self->global_best_chromosome.fitness;

	//===Compute Standard Deviation===//
	fitness_mean = compute_mean_dbl(fitness,self->num_chromosomes+1);
	fitness_std = sqrt(compute_variance_dbl(fitness, self->num_chromosomes+1));

	//===Scale===//
	for (c=0; c<self->num_chromosomes; c++){
		temp_dbl = 1.0 + (self->chromosomes[c].fitness - fitness_mean)/(2.0 * fitness_std);
		self->chromosomes[c].fitness = MAX(temp_dbl, 2*EPS);
	}
	temp_dbl = 1.0 + (self->global_best_chromosome.fitness - fitness_mean)/(2.0 * fitness_std);
	self->global_best_chromosome.fitness = MAX(temp_dbl, 2*EPS);

	return;
}

void rank_genetic_population_fitness( genetic_optimizer_t* self )
{
	unsigned int c;
	unsigned int fitness_rankings[(MAX_NUM_CHROMOSOMES+1)];
	double fitness[(MAX_NUM_CHROMOSOMES+1)];

	//===Get All Fitness Values===//
	for (c=0; c<self->num_chromosomes; c++){
		fitness[c] = self->chromosomes[c].fitness;
	}

	//===Sort Fitnesses===//
	sort_array_indices_dbl(fitness, self->num_chromosomes, fitness_rankings);

	//===Replace Fitnesses With Rankings===//
	for (c=0; c<self->num_chromosomes; c++){
		self->chromosomes[fitness_rankings[c]].fitness = c + 1;
	}
	self->global_best_chromosome.fitness = self->num_chromosomes + 1;
	
	return;
}


//should be an update global best chromosome function

void update_global_best_chromosome( genetic_optimizer_t* self )
{

	unsigned int c, best_index;
	double max_fitness;

	//===Find Best===//
	max_fitness = -(DBL_MAX-1); best_index = 0;
	for (c=0; c<self->num_chromosomes; c++){
		if (self->chromosomes[c].fitness >= max_fitness){
			max_fitness = self->chromosomes[c].fitness;
			best_index = c;
		}
	}

	//===Copy Over===//
	if ((self->chromosomes[best_index].fitness-self->global_best_chromosome.fitness)>0){
		copy_chromosome(&(self->chromosomes[best_index]), &(self->global_best_chromosome));
	}

	//===Add Fitness===//
	enqueue_ring_buffer_dbl(&(self->fitness_history), &max_fitness, 1);

	return;

}

void find_genetic_best_chromosome( genetic_optimizer_t* self,
								   chromosome_t* best_chromosome )
{
	unsigned int c, best_index;
	double max_fitness;

	//===Find Best===//
	max_fitness = -(DBL_MAX-1); best_index = 0;
	for (c=0; c<self->num_chromosomes; c++){
		if (self->chromosomes[c].fitness >= max_fitness){
			max_fitness = self->chromosomes[c].fitness;
			best_index = c;
		}
	}


	//===Copy Over===//
	if (best_chromosome != NULL){
		copy_chromosome(&(self->chromosomes[best_index]), best_chromosome);
	}
	if ((self->chromosomes[best_index].fitness-self->global_best_chromosome.fitness)>0){
		copy_chromosome(&(self->chromosomes[best_index]), &(self->global_best_chromosome));
		#if TEST_GENETIC_OPTIMIZER
			self->global_best_chromosome.calculate_fitness(&(self->global_best_chromosome));
		#endif
	}

	return;
}

void print_genetic_population_best_chromosome( genetic_optimizer_t* self )
{
	chromosome_t best;

	find_genetic_best_chromosome(self, &best);
	fprintf(stdout, "Best Chromosome:\n	");
	print_chromosome(&best, stdout);
	
	return;
}

void print_genetic_global_best_chromosome( genetic_optimizer_t* self )
{
	FILE* fout;
	fprintf(stdout, "Global Best Chromosome:\n	");
	print_chromosome(&(self->global_best_chromosome), stdout);
	fout = fopen("global_best_chromosome.dat", "w");
	print_chromosome(&(self->global_best_chromosome), fout);
	fclose(fout);

	return;
}

void select_genetic_parents_stochastic( genetic_optimizer_t* self,
										chromosome_t* parents,
										unsigned int* parent_population_indices )
{
	unsigned int c, num_parents;
	double fitness_sum, random_number, fitness_accumulated;

	//===Calculate Sum Of All Fitnesses In Entire Population===//
	fitness_sum = 0;
	for (c=0; c<self->num_chromosomes; c++){
		fitness_sum += self->chromosomes[c].fitness;
	}
	
	//===Shuffle Population===//
	shuffle_chromosomes(self->chromosomes, self->num_chromosomes, &(self->rng));
	
	//===Generate Fitness Cutoff===//
	random_number = get_random_uniform_dbl(0, (fitness_sum/((double)self->num_chromosomes)), &(self->rng));

	//===Select Parent Population Via Stochastic Univeral Sampling===//
	num_parents = 0; fitness_accumulated = 0; c = 0;
	while(num_parents < self->generation_size){

		//===Accumulate Fitness===//
		fitness_accumulated += self->chromosomes[c].fitness;

		//===Add Parent To Population A Number Of Times According To Fitness===//
		while(fitness_accumulated > random_number){
			parent_population_indices[num_parents] = c;
			random_number += fitness_sum/((double)self->num_chromosomes);
			num_parents++;
			if (num_parents >= self->generation_size) break;
		}

		//===Increment Parent===//
		c++;
	}

	//===Copy Over Parents===//
	for (c=0; c<num_parents; c++){
		copy_chromosome(&(self->chromosomes[parent_population_indices[c]]), &(parents[c]));
	}

	return;
}

void find_unique_elites( chromosome_t* chromosomes,
						 unsigned int num_chromosomes,
						 chromosome_t* elites,
						 unsigned int num_elites )
{
	unsigned int c, e, found_elites, equal;

	if (num_elites == 0){
		fprintf(stderr, "Error:: Num Elites Is Zero! In Function -- find_unique_elites!\n");
		quit();
	}

	//===Sort Parents By Fitness===//
	sort_chromosomes_by_fitness(chromosomes, num_chromosomes);

	//===Copy 1st Elite===//
	copy_chromosome(&(chromosomes[0]), &(elites[0]));
	if (num_elites == 1) return;

	//===Find All Other Elites===//
	found_elites = 1; c = 1;
	while(found_elites < num_elites){
		
		//===Is This Chromosome Equal To Any Of The Elites===//
		equal = 1;
		for (e=0; e<found_elites; e++){
			equal = equal && chromosomes_are_equal(&(elites[e]), &(chromosomes[c]));
		}
		if (!equal){
			copy_chromosome(&(chromosomes[c]), &(elites[found_elites]));
			found_elites++;
		}
		
		//===Increment Chromosome===//
		c++;
	}

	return;
}
							
void breed_genetic_population( genetic_optimizer_t* self )
{
	unsigned int c, parent, num_children;
	unsigned int parent_population_indices[(MAX_NUM_CHROMOSOMES+1)];
	#if ALLOW_GLOBAL_BEST_REPLACEMENT
		unsigned int equal, least_fit_chromosome_index;
	#endif
		chromosome_t parents[(MAX_NUM_CHROMOSOMES+1)], children[(MAX_NUM_CHROMOSOMES+1)];
	#if ALLOW_ELITES
		chromosome_t elite_population[(MAX_NUM_CHROMOSOMES+1)];
	#endif
	chromosome_t parent_population[(MAX_NUM_CHROMOSOMES+1)];
	chromosome_t child_population[(MAX_NUM_CHROMOSOMES+1)];

	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Genetic Optimizer Is NULL! In Function -- breed_genetic_population_stud_2!\n");
		quit();
	}

	//===Select Parents Via Stochastic Universal Sampling===//
	select_genetic_parents_stochastic(self, parent_population, parent_population_indices);

	//===Run Elite Selection===//
	num_children = 0;
	#if ALLOW_ELITES	
		//===Find Unique Elites===//
		if (self->num_elites == 1){
			find_most_fit_chromosome(parent_population, self->generation_size, &(child_population[0]), NULL);
			num_children++;
		}
		else if (self->num_elites > 1){

			find_unique_elites(parent_population, self->generation_size, elite_population, self->num_elites);

			//===Copy Elites To Children===//
			while(num_children < self->num_elites){

				//===Copy Elite===//
				copy_chromosome(&(elite_population[num_children]), &(child_population[num_children]));		

				//===Increment Children===//
				num_children++;
			}		

		}
		else{
			num_children = 0;
		}
	#endif

	#if ALLOW_GLOBAL_BEST_REPLACEMENT
		//===Check If Global Best Is The Same As Any In The Population===//
		equal = 0;
		for (c=0; c<self->generation_size; c++){
			equal = chromosomes_are_equal(&(self->global_best_chromosome), &(parent_population[c]));
		}
		//===If Not Replace Least Fit Parent With Global Best===//
		if (equal == 0){
			least_fit_chromosome_index = UINT_MAX;
			find_least_fit_chromosome(parent_population, self->generation_size, NULL, &least_fit_chromosome_index);
			if (least_fit_chromosome_index != UINT_MAX){
				copy_chromosome(&(self->global_best_chromosome), 
								&(parent_population[least_fit_chromosome_index]));			
			}
		}
	#endif

	#if ALLOW_STUD	

		//===Sort Parents By Fitness===//
		sort_chromosomes_by_fitness(parent_population, self->generation_size);

		//===Mate The Stud===//
		parent = 1;
		copy_chromosome(&(parent_population[0]), &(parents[0]));		
		while( num_children < self->generation_size ){

			#if ALLOW_MULTIPLE_PARENTS
				//===Select Other Parents===//
				for (c=1; c<self->num_parents; c++){
					parent = get_random_uniform(1, self->generation_size, &(self->rng));
					copy_chromosome(&(parent_population[parent]), &(parents[c]));
				}
			#else
				//===Select 2nd Parent===//
				copy_chromosome(&(parent_population[parent]), &(parents[1]));
				if (chromosomes_are_equal(&(parents[0]), &(parents[1]))){	
					//===Increment Parent===//
					parent = ((parent + 1) % (self->generation_size));
					continue;
				}		
			#endif

			//===Breed Chromosomes===//
			breed_chromosomes_fuzzy(parents, children, self->num_parents, &(self->rng));
			
			//===Save And Increment Children===//
			if (num_children < self->num_chromosomes){
				copy_chromosome(&(children[0]), &(child_population[num_children++]));
			}
			else{
				break;
			}

			//===Increment Parent===//
			parent = ((parent + 1) % (self->generation_size));
		}

		//===Mutate Children That Are Not Elites===//
		#if ALLOW_ELITES
			for (c=self->num_elites; c<self->generation_size; c++){
				mutate_chromosome_genes(&(child_population[c]), &(self->rng));
			}
		#else
			for (c=0; c<self->generation_size; c++){
				mutate_chromosome_genes(&(child_population[c]), &(self->rng));
			}
		#endif

	#else

		//===Mate 2 by 2===//
		parent = 0;
		while( num_children < self->generation_size ){

			#if ALLOW_MULTIPLE_PARENTS
				//===Select Parents===//
				for (c=0; c<self->num_parents; c++){
					parent = get_random_uniform(0, self->generation_size, &(self->rng));
					copy_chromosome(&(parent_population[parent]), &(parents[c]));
				}
			#else
				//===Select 1st Parent===//
				copy_chromosome(&(parent_population[parent]), &(parents[0]));
				parent = ((parent + 1) % (self->generation_size));

				//===Select 2nd Parent===//
				copy_chromosome(&(parent_population[parent]), &(parents[1]));
				parent = ((parent + 1) % (self->generation_size));
			#endif

			//===Breed Chromosomes===//
			breed_chromosomes_fuzzy(parents, children, self->num_parents, &(self->rng));

			//===Save And Increment Children===//
			if (num_children < self->num_chromosomes){
				copy_chromosome(&(children[0]), &(child_population[num_children++]));
			}
			else{
				break;
			}
			if (num_children < self->num_chromosomes){
				copy_chromosome(&(children[1]), &(child_population[num_children++]));
			}
			else{
				break;
			}

		}
		//===Mutate Children That Are Not Elites===//
		#if ALLOW_ELITES
			for (c=self->num_elites; c<self->generation_size; c++){
				mutate_chromosome_genes(&(child_population[c]), &(self->rng));
			}
		#else
			for (c=0; c<self->generation_size; c++){
				mutate_chromosome_genes(&(child_population[c]), &(self->rng));
			}
		#endif

	#endif

	//===Sort Population By Fitness===//
	sort_chromosomes_by_fitness(self->chromosomes, self->num_chromosomes);

	//===Copy Children To Least Fit Parents===//
	for (c=0; c<self->generation_size; c++){
		copy_chromosome(&(child_population[c]), &(self->chromosomes[self->num_chromosomes-c-1]));
	}



	return;
}

void calculate_genetic_population_fitness( genetic_optimizer_t* self )
{
	unsigned int c;

	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Genetic Optimizer Is NULL! In Function -- calculate_genetic_population_fitness!\n");
		quit();
	}

	//===Calculate Fitness For Each Chromosome===//
	for (c=0; c<self->num_chromosomes; c++){
		self->chromosomes[c].calculate_fitness(&(self->chromosomes[c]));
	}
	self->global_best_chromosome.calculate_fitness(&(self->global_best_chromosome));

	//===Sort Population By Fitness===//
	sort_chromosomes_by_fitness(self->chromosomes, self->num_chromosomes);

	
	return;
}

//================================================================================================//
//========================================TEST METHODS============================================//
//================================================================================================//

void evaluate_sphere_benchmark_fitness( chromosome_t* self )
{
	unsigned int g;
	double fitness;
	double temp;

	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Chromosome Is NULL! In Function -- evaluate_sphere_benchmark_fitness!\n");
		quit();
	}

	fitness = 0;
	for (g=0; g<self->num_genes; g++){
		temp = self->genes[g].allele - 5.0;
		fitness += (temp*temp); 
	}
	fitness *= -1.0;
	self->fitness = fitness;

	//===Check For Negative Zero===//
	if (fabs(self->fitness) < 2*EPS){
		self->fitness = 0.0;
	}

	return;
}

void evaluate_ackley_benchmark_fitness( chromosome_t* self )
{

	unsigned int g;
	double temp1, temp2;
	double fitness;

	//===Check For Errors===//
	if (self == NULL){
		fprintf(stderr, "Error:: Chromosome Is NULL! In Function -- evaluate_ackley_benchmark_fitness!\n");
		quit();
	}

	//===Calcuate Sums===//
	temp1 = 0;
	for (g=0; g<self->num_genes; g++){
		temp1 += (self->genes[g].allele*self->genes[g].allele);
	}
	temp1 /= ((double)(self->num_genes));
	temp2 = 0;
	for (g=0; g<self->num_genes; g++){
		temp2 += cos(2.0*M_PI*self->genes[g].allele);
	}
	temp2 /= ((double)(self->num_genes));

	//===Make Sum===//
	fitness = 0;
	for (g=0; g<self->num_genes; g++){
		fitness += (20.0 + exp(1) - 20.0 * exp(-0.2*temp1) - exp(temp2));
	}
	fitness *= -1.0;
	self->fitness = fitness;

	//===Check For Negative Zero===//
	if (fabs(self->fitness) < 2*EPS){
		self->fitness = 0.0;
	}

	return;
}

void evaluate_griewank_benchmark_fitness( chromosome_t* self )
{

	unsigned int g;
	double temp1, temp2;
	double fitness;

	//===Calcuate Sum===//
	temp1 = 0;
	for (g=0; g<self->num_genes; g++){
		temp1 += (self->genes[g].allele*self->genes[g].allele);
	}
	temp1 /= ((4000.0));

	//===Calculate Product===//
	temp2 = 1;
	for (g=0; g<self->num_genes; g++){
		temp2 *= cos((self->genes[g].allele/sqrt(((double)(g+1)))));
	}

	//===Make Sum===//
	fitness = 0;
	for (g=0; g<self->num_genes; g++){
		fitness += (1.0 + temp1 - temp2);
	}
	fitness *= -1.0;
	self->fitness = fitness;

	//===Check For Negative Zero===//
	if (fabs(self->fitness) < 2*EPS){
		self->fitness = 0.0;
	}

	return;
}

void evaluate_schwefel_sine_benchmark_fitness( chromosome_t* self )
{

	unsigned int g;
	double fitness;

	//===Make Sum===//
	fitness = 0;
	for (g=0; g<self->num_genes; g++){
		fitness += (-1.0 * self->genes[g].allele * sin(sqrt(fabs(self->genes[g].allele))));
	}
	self->fitness = fabs(fitness + 12965.5);

	//===Check For Negative Zero===//
	if (fabs(self->fitness) < 2*EPS){
		self->fitness = 0.0;
	}

	return;
}

void evaluate_simple_quadratic_benchmark_fitness( chromosome_t* self )
{
	unsigned int g;
	double fitness;

	//===Make Sum===//
	fitness = 0;
	for (g=0; g<self->num_genes; g++){
		fitness += (self->genes[g].allele * self->genes[g].allele);
	}
	fitness -= 1.0; 
	fitness *= -1.0;
	self->fitness = fitness;

	//===Check For Negative Zero===//
	if (fabs(self->fitness) < 2*EPS){
		self->fitness = 0.0;
	}
	return;
}

void evaluate_fit_cosine_benchmark_fitness( chromosome_t* self )
{
	double fitness;

	//===Make Cosine===//
	fitness = (2.0 - self->genes[0].allele) * cos(2.0*M_PI*(100.0 - self->genes[1].allele)/1600.0  + ((M_PI/2.0) - self->genes[2].allele) );
	fitness *= -1.0;
	self->fitness = fitness;

	//===Check For Negative Zero===//
	if (fabs(self->fitness) < 2*EPS){
		self->fitness = 0.0;
	}
	return;
}

void test_genetic_optimizer()
{
	unsigned int iteration, keep_running, iteration_limit;
	unsigned int g, c, num_genes, num_parents, num_chromosomes, num_elites, generation_size;
	double fitness_mean, fitness_variance, fitness_mean_threshold, fitness_variance_threshold;
	double gene_bounds[2*MAX_NUM_GENES];
	genetic_optimizer_t* optimizer;

	//===Mallocs===//
	optimizer = malloc(sizeof(genetic_optimizer_t));

	//===Initialize Optimizer===//
	num_genes = 5;
	num_chromosomes = 20;
	num_elites = 0;
	num_parents = 3;
	generation_size = num_chromosomes;
	for (g=0; g<num_genes; g++){
		gene_bounds[2*g + 0] = -500.0;
		gene_bounds[2*g + 1] = 500.0;
	}

	//===Initialize Optimizer===//
	initialize_genetic_optimizer(optimizer,gene_bounds, NULL, NULL, num_genes, num_chromosomes,num_elites,num_parents,generation_size);

	//===Set Fitness Function===// -- should do this in the initializer
	for (c=0; c<num_chromosomes; c++){
		//optimizer->chromosomes[c].calculate_fitness = evaluate_sphere_benchmark_fitness;
		//optimizer->chromosomes[c].calculate_fitness = evaluate_ackley_benchmark_fitness;
		//optimizer->chromosomes[c].calculate_fitness = evaluate_griewank_benchmark_fitness;
		optimizer->chromosomes[c].calculate_fitness = evaluate_schwefel_sine_benchmark_fitness;
		//optimizer->chromosomes[c].calculate_fitness = evaluate_simple_quadratic_benchmark_fitness;
		//optimizer->chromosomes[c].calculate_fitness = evaluate_fit_cosine_benchmark_fitness;

	}
	//optimizer->global_best_chromosome.calculate_fitness = evaluate_sphere_benchmark_fitness;
	//optimizer->global_best_chromosome.calculate_fitness = evaluate_ackley_benchmark_fitness;
	//optimizer->global_best_chromosome.calculate_fitness = evaluate_griewank_benchmark_fitness;
	optimizer->global_best_chromosome.calculate_fitness = evaluate_schwefel_sine_benchmark_fitness;
	//optimizer->global_best_chromosome.calculate_fitness = evaluate_simple_quadratic_benchmark_fitness;
	//optimizer->global_best_chromosome.calculate_fitness = evaluate_fit_cosine_benchmark_fitness;

	//===Run Optimizer===//
	fitness_mean_threshold = 2.0*EPS; fitness_variance_threshold = 10.0*EPS;
	keep_running = 1; iteration = 0; iteration_limit = 1000;
	while(keep_running){

		//===Calculate Fitness===//
		fprintf(stdout, "Iteration: %d\n", iteration);
		calculate_genetic_population_fitness(optimizer);

		//===Update Best===//
		update_global_best_chromosome(optimizer);

		//===Print Best===//
		print_genetic_global_best_chromosome(optimizer);

		//===Print===//
		print_genetic_population(optimizer, stdout);
		newline();

		//===Scale Fitness===//
		scale_genetic_population_fitness(optimizer);

		//===Rank Population===//
		if (iteration < iteration_limit/10){
			rank_genetic_population_fitness(optimizer);
		}

		//===Breed===//	
		breed_genetic_population(optimizer);

		//===Incremement===//
		iteration++;

		//if (iteration > 1) quit();

		//===Test Break===//
		if (iteration > optimizer->fitness_history.queue_length){

			//===Calculate Mean Of Best Fitnesses===//
			fitness_mean = calculate_ring_buffer_mean_dbl(&(optimizer->fitness_history));
			//fprintf(stdout, "Mean: %lf\n", fitness_mean);

			//===Calculate Variances Of Best Fitnesses===//
			fitness_variance = calculate_ring_buffer_variance_dbl(&(optimizer->fitness_history));

			//===Compare===//
			if (fabs(optimizer->global_best_chromosome.fitness - fitness_mean) < fitness_mean_threshold){
				if (fabs(fitness_variance) < fitness_variance_threshold){
					keep_running = 0;
				}
			}
			
		}	

		//===Hard Iteration Limit===//
		if (iteration >= iteration_limit){
			keep_running = 0;
		}

	}	

	//===Clean Up===//
	free(optimizer);

	return;
}
