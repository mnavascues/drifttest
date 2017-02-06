/*
 Copyright INRA
 author: Miguel Navascués (2017)

 based on code by Renaud Vitalis

 This file is part of DriftTest.

 DriftTest is a computer program whose purpose is to is to detect
 loci under selection based on changes in allele frequencies
 between two time samples of the same population.

 DriftTest is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 DriftTest is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <errno.h>           // system error numbers
#include <stdio.h>           // standrad input/output
#include <stdlib.h>          // general purpose standard library
#include <getopt.h>          // for parsing command line arguments
#include <math.h>            // for computing common mathematical operations
#include <gsl/gsl_rng.h>     // GSL: random number generation
#include <gsl/gsl_randist.h> // GSL: random number distributions

#include "defs.h"
#include "popgen.h"

int main (int argc, char *argv[])
{
  printf ("\ndrifttest\n\n");

  //declare random number generator and set type to “Mersenne Twister” (MT19937)
  gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
  seed = 0;

  unsigned int Ne;
  double Fis, obs_Fst;
  double p;
  unsigned int genotype_counts[2][3];

  //declare variables for parsing command line arguments
  int opt = 0;
  int long_index = 0;
  static const struct option long_options[] =
  {
    {"seed",required_argument,NULL,1},
    {"tau",required_argument,NULL,2},
    {"maf",required_argument,NULL,3},
    {"help",no_argument,NULL,201},
    {"version",no_argument,NULL,202},
    {"test",no_argument,NULL,203},
    {NULL,0,NULL,0}
  };

  // parse command line arguments (with getopt_long_only function from getopt.h)
  program_name = argv[0];
  while ((opt = getopt_long_only(argc, argv,"",long_options,&long_index)) != -1)
  {
    switch (opt) {
      case 1 :
        seed = atoi(optarg);
        printf ("Seed for random number generator: ’%lu’ \n", seed);
        break;
      case 2 :
        tau = atoi(optarg);
        printf ("Number of generations between samples (tau): ’%d’ \n", tau);
        break;
      case 3 :
        maf = atof(optarg);
        printf ("Minimum allele frequency (maf): ’%f’ \n", maf);
        break;
      case 201 :
        print_usage();
        break;
      case 202 :
        print_version();
        exit(EXIT_SUCCESS);
      case 203 :
	seed = 123456;
	tau  = 10;
	maf  = 0.1;
	gsl_rng_set(r,seed);
        print_test(r);
        exit(EXIT_SUCCESS);
      default :
        print_usage();
        exit(EXIT_FAILURE);
    }
  }
  if (optind < argc)
  {
    fprintf(stderr, "Error! Unrecognized option(s) ");
    while (optind < argc)
      printf ("`%s' ", argv[optind++]);
    printf ("\n");
    print_usage();
    exit(EXIT_FAILURE);
  }
  if (seed == 0) {
    fprintf(stderr, "Warning: seed has value ’%lu’; using default GSL seed\n", seed);
  }
  if (tau <= 0) {
    fprintf(stderr, "tau: '%d' \n",tau);
    fprintf(stderr, "Error! The value of option -tau has to be positive\n");
    exit(EXIT_FAILURE);
  }

  // set seed for random number generator
  //printf ("Random number generator: ’%s’ \n", gsl_rng_name (r));
  gsl_rng_set(r,seed);













  Ne = 111;
  Fis = 0.8;
  genotype_counts[0][0] = 20;
  genotype_counts[0][1] = 15;
  genotype_counts[0][2] = 65;
  genotype_counts[1][0] = 95;
  genotype_counts[1][1] = 1;
  genotype_counts[1][2] = 4;


  printf ("Ne: ’%d’ \n", Ne);
  printf ("\n");
  printf ("Fis: ’%f’ \n", Fis);
  printf ("\n");
  printf ("start AA: ’%d’ \n", genotype_counts[0][0]);
  printf ("start AB: ’%d’ \n", genotype_counts[0][1]);
  printf ("start BB: ’%d’ \n", genotype_counts[0][2]);
  printf ("\n");
  printf ("end AA: ’%d’ \n", genotype_counts[1][0]);
  printf ("end AB: ’%d’ \n", genotype_counts[1][1]);
  printf ("end BB: ’%d’ \n", genotype_counts[1][2]);



  obs_Fst = FST_2pop_from_genotypes_counts(genotype_counts);
  printf ("\n");
  printf ("Fst: ’%f’ \n", obs_Fst);
  nbr_simuls = NBR_SIMULS;
  p = p_value(r, Ne, Fis, nbr_simuls, genotype_counts);
  printf ("\n");
  printf ("p: ’%f’ \n", p);

 








  



  gsl_rng_free (r);
  printf ("\ndrifttest has finished. Good bye!\n\n");
  return (EXIT_SUCCESS);
}





// PRINT USAGE 

void print_usage()
{
  printf("usage: %s [ options ]\n",program_name);
  printf("valid options are :\n");
  printf("-help\t\t\t print this message\n");
  printf("-version\t\t print version\n");
  printf("-test\t\t\t print results of a set of tests on the code\n");
  printf("-seed\t\t\t initial seed for the random number generator (default: GSL default)\n");
  exit(EXIT_SUCCESS);
}

// PRINT VERSION 

void print_version()
{
  printf("You are running version %s\n",VERSION);
}

// PRINT TESTS 

void print_test(const gsl_rng * r)
{

  printf("You are running a set of tests to verify the integrity of drifttest code\n\n");

  if (seed != 123456)
  {
    seed = 123456; 
    printf("Seed was set to %lu\n\n", seed);
  }

   
  // drift_sim
  double allele_freq;
  allele_freq = drift_sim(r, 0.5, 10000);
  printf ("Initial allele frequency: 0.5\nTime of drift: %d generations\nStrength of drift: Ne=10000\nFinal allele frequency: %f\n\n", tau, allele_freq);
  if (allele_freq!= 0.4852) fprintf(stderr, "Warning: final frequency should be 0.4852 for seed=123456\n");
  allele_freq = drift_sim(r, 0.5, 10);
  printf ("Initial allele frequency: 0.5\nTime of drift: %d generations\nStrength of drift: Ne=10\nFinal allele frequency: %f\n\n", tau, allele_freq);
  if (allele_freq!= 0.2) fprintf(stderr, "Warning: final frequency should be 0.2 for seed=123456\n");



}