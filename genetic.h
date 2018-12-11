/** @file genetic.h
*   @brief Contains functions for genetically optimizing parameters.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date August 2016
*  @bug No known bugs
*/

#ifndef GENETIC_H
#define GENETIC_H

//================================================================================================//
//===================================STANDARD INCLUDES============================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"
#include "ring.h"
#include "random_numbers.h"



//================================================================================================//
//========================================MACROS==================================================//
//================================================================================================//

#define ALLOW_DOMINANT_GENETICS 0
#define ALLOW_STUD 1
#define ALLOW_ELITES 0
#define ALLOW_MUTATION_RATE_RECOMBINATION 1
#define ALLOW_MUTATION_MAGNITUDE_RECOMBINATION 1
#define ALLOW_MULTIPLE_PARENTS 1
#define ALLOW_GAUSSIAN_MUTATION 0
#define ALLOW_DOMAIN_CENTERED_MUTATION 0
#define TEST_GENETIC_OPTIMIZER 0
#define ALLOW_GLOBAL_BEST_REPLACEMENT 1


//================================================================================================//
//======================================DATA STRUCTURES===========================================//
//================================================================================================//

//================================================================================================//
/** @struct gene_t
*   @brief This structure is the typedef for the gene_t object.
*/
//================================================================================================//
typedef struct gene_s gene_t;
typedef struct gene_s{
	double mutation_rate;
	double mutation_magnitude;
	double bounds[2];
	double allele;
} gene_t;


//================================================================================================//
/** @struct chromosome_t
*   @brief This structure is the typedef for the chromosome_t object.
*/
//================================================================================================//
typedef struct chromosome_s chromosome_t;
typedef struct chromosome_s{
	unsigned int num_genes;
	unsigned int dominant;
	double fitness;
	gene_t genes[MAX_NUM_GENES];
	void (*calculate_fitness)(chromosome_t*);
} chromosome_t;


//================================================================================================//
/** @struct genetic_optimizer_t
*   @brief This structure is the typedef for the genetic_optimizer_t object.
*/
//================================================================================================//
typedef struct genetic_optimizer_s genetic_optimizer_t;
typedef struct genetic_optimizer_s{
	unsigned int num_chromosomes;
	unsigned int num_elites;
	unsigned int num_parents;
	unsigned int generation_size;
	ring_buffer_t fitness_history;
	chromosome_t global_best_chromosome;
	chromosome_t chromosomes[MAX_NUM_CHROMOSOMES];
	random_number_generator_t rng;
	void (*breed)(genetic_optimizer_t*);
} genetic_optimizer_t;



//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//================================================================================================//
/**
* @brief This function initializes a gene_t object.
*
* @param[in,out] gene_t* self
* @param[in] double* bounds
* @param[in] double allele
* @param[in] double mutation_rate
* @param[in] double mutation_magnitude
*
* @return NONE
*/
//================================================================================================//
void initialize_gene(gene_t*,double*,double,double,double);


//================================================================================================//
/**
* @brief This function copies a gene_t object.
*
* @param[in] gene_t* original
* @param[out] gene_t* copy
*
* @return NONE
*/
//================================================================================================//
void copy_gene(gene_t*,gene_t*);


//================================================================================================//
/**
* @brief This function randomly chooses the gene_t object's attributes.
*
* NOTE: How it chooses is dependent upon the user chosen flags.
*
* @param[in,out] gene_t* original
* @param[in] random_number_generator_t* rng
*
* @return NONE
*/
//================================================================================================//
void mutate_gene(gene_t*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function initializes a chromosome_t object.
*
* @param[in,out] chromosome_t* self
* @param[in] double* gene_bounds
* @param[in] double* mutation_magnitudes
* @param[in] double* mutation_rates
* @param[in] unsigned int num_genes
* @param[in] random_number_generator_t* rng
*
* @return NONE
*/
//================================================================================================//
void initialize_chromosome(chromosome_t*,double*,double*,double*,unsigned int,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function prints a chromosome to a file.
*
* @param[in] chromosome_t* self
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_chromosome(chromosome_t*,FILE*);


//================================================================================================//
/**
* @brief This function copies a chromosome_t object.
*
* @param[in] chromosome_t* original
* @param[out] chromosome_t* copy
*
* @return NONE
*/
//================================================================================================//
void copy_chromosome(chromosome_t*,chromosome_t*);


//================================================================================================//
/**
* @brief This function randomly chooses all the chromosomes's gene's attributes.
*
* NOTE: How it chooses is dependent upon the user chosen flags.
*
* @param[in,out] chromosome_t* original
* @param[in] random_number_generator_t* rng
*
* @return NONE
*/
//================================================================================================//
void mutate_chromosome_genes(chromosome_t*,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function tests if two chromosome objects have all the same fields.
*
* @param[in] chromosome_t* chromosome1
* @param[in] chromosome_t* chromosome2
*
* @return unsigned int are_equal
*/
//================================================================================================//
unsigned int are_chromosomes_equal(chromosome_t*,chromosome_t*);


//================================================================================================//
/**
* @brief This function breeds two children from a number of parent chromosomes.
*
* NOTE: This function uses a triangular distribution to linearly combine genes.
*
* @param[in] chromosome_t* parents
* @param[out] chromosome_t* children
* @param[in] unsigned int num_parents
* @param[in] random_number_generator_t* rng
*
* @return NONE
*/
//================================================================================================//
void breed_chromosomes_fuzzy(chromosome_t*,chromosome_t*,unsigned int,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function sorts a group of chromosomes in decreasing order according to fitnesses.
*
* @param[in,out] chromosome_t* chromosomes
* @param[in] unsigned int num_chromosomes
*
* @return NONE
*/
//================================================================================================//
void sort_chromosomes_by_fitness(chromosome_t*,unsigned int);


//================================================================================================//
/**
* @brief This function shuffles an array of chromosomes.
*
* @param[in,out] chromosome_t* chromosomes
* @param[in] unsigned int num_chromosomes
* @param[in] random_number_generator_t* rng
*
* @return NONE
*/
//================================================================================================//
void shuffle_chromosomes(chromosome_t*,unsigned int,random_number_generator_t*);


//================================================================================================//
/**
* @brief This function initializes a genetic_optimizer_t object.
*
* @param[in,out] genetic_optimizer_t* self
* @param[in] double* gene_bonds
* @param[in] double* mutation_magntiudes
* @param[in] double* mutation_rates
* @param[in] unsigned int num_genes
* @param[in] unsigned int num_chromosomes
* @param[in] unsigned int num_elites
* @param[in] unsigned int num_parents
* @param[in] unsigned int generation_size
*
* @return NONE
*/
//================================================================================================//
void initialize_genetic_optimizer(genetic_optimizer_t*,double*,double*,double*,unsigned int,
								   unsigned int,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function prints each chromosome in the population to a file.
*
* @param[in] genetic_optimizer_t* self
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_genetic_population(genetic_optimizer_t*,FILE*);


//================================================================================================//
/**
* @brief This function sigma scales the fitnesses of each chromosome in the population.
*
* @param[in,out] genetic_optimizer_t* self
*
* @return NONE
*/
//================================================================================================//
void scale_genetic_population_fitness(genetic_optimizer_t*);


//================================================================================================//
/**
* @brief This function replaces fitness of each population chromosome with their relative rank.
*
* @param[in,out] genetic_optimizer_t* self
*
* @return NONE
*/
//================================================================================================//
void rank_genetic_population_fitness(genetic_optimizer_t*);


//================================================================================================//
/**
* @brief This function updates the best fit chromosome in the population and adds fitness to buffer.
*
* @param[in] genetic_optimizer_t* self
*
* @return NONE
*/
//================================================================================================//
void update_global_best_chromosome(genetic_optimizer_t*);


//================================================================================================//
/**
* @brief This function prints the best fit chromosome in the population.
*
* @param[in] genetic_optimizer_t* self
*
* @return NONE
*/
//================================================================================================//
void print_genetic_population_best_chromosome(genetic_optimizer_t*);


//================================================================================================//
/**
* @brief This function prints the best fit chromosome found so far.
*
* @param[in] genetic_optimizer_t* self
*
* @return NONE
*/
//================================================================================================//
void print_genetic_global_best_chromosome(genetic_optimizer_t*);


//================================================================================================//
/**
* @brief This function selects a population of parents via stochastic universal sampling.
*
* @param[in] genetic_optimizer_t* self
* @param[in,out] chromosome_t* parents
* @param[in,out] unsigned int* parent_population_indices
*
* @return NONE
*/
//================================================================================================//
void select_genetic_parents_stochastic(genetic_optimizer_t*,chromosome_t*,unsigned int*);


//================================================================================================//
/**
* @brief This function selects the top chromosomes as elites.
*
* @param[in] chromosome_t* chromosomes
* @param[in] unsigned int num_chromosomes
* @param[out] chromosome_t* elites
* @param[in] unsigned int num_elites
*
* @return NONE
*/
//================================================================================================//
void find_unique_elites(chromosome_t*,unsigned int,chromosome_t*,unsigned int);


//================================================================================================//
/**
* @brief This function creates the next generation of chromosomes.
*
* @param[in,out] genetic_optimizer_t* self
*
* @return NONE
*/
//================================================================================================//
void breed_genetic_population(genetic_optimizer_t*);


//================================================================================================//
/**
* @brief This function tests the genetic optimizer against some benchmark functions.
*
* @return NONE
*/
//================================================================================================//
void test_genetic_optimizer();


#endif //GENETIC_H//
