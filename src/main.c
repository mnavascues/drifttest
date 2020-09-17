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
#include <string.h>
#include <getopt.h>          // for parsing command line arguments
#include <math.h>            // for computing common mathematical operations
#include <gsl/gsl_rng.h>     // GSL: random number generation
#include <gsl/gsl_randist.h> // GSL: random number distributions
#include <omp.h>             // OpenMP (parallel computing)


#include "main.h"
#include "popgen.h"
#include "test.h"
#include "readdata.h"
#include "writeresults.h"

int main (int argc, char *argv[])
{
  printf ("\n\n\nThis is DriftTest by Miguel Navascués\n");
  printf("(version %s)\n\n\n",VERSION);
  printf ("DriftTest is a computer program whose purpose is to is to detect\n");
  printf ("loci under selection based on changes in allele frequencies\n");
  printf ("between two time samples of the same population.\n\n\n");

  unsigned long seed  = 0;
  char data_filename[106];
    strcpy(data_filename,"data/data.txt");
  char results_filename[124];
    strcpy(results_filename,"results/locus_by_locus");
  char multilocus_results_filename[120];
    strcpy(multilocus_results_filename,"results/multilocus");
  int n_threads = 0;
  int locus;
  int user_supplied_fis = 0;
  data_struct data;
  global_result_struct global_result;
  unsigned int one_locus_genotype_counts[2][3];
  gsl_rng * r;

  //declare variables for parsing command line arguments
  int opt = 0;
  int long_index = 0;
  static const struct option long_options[] =
  {
    {"seed",required_argument,NULL,1},
    {"tau",required_argument,NULL,2},
    {"maf",required_argument,NULL,3},
    {"infile",required_argument,NULL,4},
    {"outfile",required_argument,NULL,5},
    {"fis",required_argument,NULL,6},
    {"threads",required_argument,NULL,7},
    {"help",no_argument,NULL,201},
    {"version",no_argument,NULL,202},
    {"test",no_argument,NULL,203},
    {NULL,0,NULL,0}
  };

  //default values
  maf = 0.0;
  fis = 0.0;

  // parse command line arguments (with getopt_long_only function from getopt.h)
  program_name = argv[0];
  while ((opt = getopt_long_only(argc, argv,"",long_options,&long_index)) != -1)
  {
    switch (opt) {
      case 1 :
        seed = atoi(optarg);
        printf ("Seed for random number generator: %lu\n", seed);
        break;
      case 2 :
        tau = atoi(optarg);
        printf ("Number of generations between samples (tau): %d\n", tau);
        break;
      case 3 :
        maf = atof(optarg);
        printf ("Minimum allele frequency threshold (maf): %f\n", maf);
        break;
      case 4 :
        strcpy(data_filename,optarg);
        printf ("Input data file: %s\n", data_filename);
        break;
      case 5 :
        strcpy(results_filename,optarg);
        strcpy(multilocus_results_filename,optarg);
        strcat(results_filename,"locus_by_locus");
        strcat(multilocus_results_filename,"multilocus");
        printf ("Output results file: %s\n", results_filename);
        printf ("Multilocus output results file: %s\n", multilocus_results_filename);
        break;
      case 6 :
        fis = atof(optarg);
        user_supplied_fis = 1;
        printf ("Inbreeding coefficient (Fis): %f\n", fis);
        break;
      case 7 :
        n_threads = atoi(optarg);
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
        r = gsl_rng_alloc (gsl_rng_mt19937);
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
    fprintf(stderr, RED "Error! Unrecognized option(s) " RESET);
    while (optind < argc)
      printf ("`%s' ", argv[optind++]);
    printf ("\n");
    print_usage();
    exit(EXIT_FAILURE);
  }
  if (seed == 0) {
    fprintf(stderr, RED "Warning: seed has value ’%lu’; using default GSL seed\n" RESET, seed);
  }
  if (tau <= 0) {
    fprintf(stderr, RED "tau: '%d' \n" RESET,tau);
    fprintf(stderr, RED "Error! The value of option -tau has to be positive\n" RESET);
    exit(EXIT_FAILURE);
  }
  if (maf < 0 || maf > 0.5) {
    fprintf(stderr, RED "tau: '%f' \n" RESET,maf);
    fprintf(stderr, RED "Error! The value of option -maf has to be a number between 0 and 0.5\n" RESET);
    exit(EXIT_FAILURE);
  }
  if (fis < -1 || fis > 1) {
    fprintf(stderr, RED "fis: '%f' \n",fis);
    fprintf(stderr, RED "Warning! The value of option -fis has to be a number between 0 and 1\n" RESET);
    fprintf(stderr, RED "Reverting to default behaviour, Fis estimated from data\n" RESET);
    user_supplied_fis = 0;
  } else if (fis < 0.0) {
    fprintf(stderr, RED "fis: '%f' \n",fis);
    fprintf(stderr, RED "Warning! The value of option -fis has to be a number between 0 and 1\n" RESET);
    fprintf(stderr, RED "A value of Fis=0.0 (random mating) will be used instead of a negative Fis\n" RESET);
    fis=0.0;
    user_supplied_fis = 1;
  }
  if (n_threads < 0) {
    printf(RED "Error! The value of option -threads has to be positive\n" RESET);
    exit(EXIT_FAILURE);
  }
  if (n_threads > omp_get_max_threads()) {
    printf(RED "Error! The value of option -threads has to be less than the maximum number of threads available\n" RESET);
    exit(EXIT_FAILURE);
  }
  if (n_threads == 0) {
    n_threads = omp_get_max_threads();  // omp_get_max_threads() returns the same value whether executing from a serial or parallel region
  }
  omp_set_num_threads(n_threads);


  read_data(&data,data_filename);

  //F_statistics(&data,&global_result);
  //printf ("Fst = %f\n",global_result.Fst);
  //printf ("Fis = %f\n",global_result.Fis);
  //printf ("Ne  = %d\n",global_result.Ne);

  F_statistics_bootstrap(&data,&global_result);
  printf ("Fst   = %f, 95CI(%f,%f)\n",global_result.Fst,global_result.Fst_lower,global_result.Fst_upper);
  printf ("Fis   = %f, 95CI(%f,%f)\n",global_result.Fis,global_result.Fis_lower,global_result.Fis_upper);
  printf ("Ne    = %d, 95CI(%d,%d)\n",global_result.Ne,global_result.Ne_lower,global_result.Ne_upper);
  printf ("sigma = %f, 95CI(%f,%f)\n",global_result.sigma,global_result.sigma_lower,global_result.sigma_upper);



  global_result.FST = (double *) malloc(data.nbr_loci * sizeof(double));
  global_result.pvalue = (double *) malloc(data.nbr_loci * sizeof(double));

#pragma omp parallel for schedule(guided) private(one_locus_genotype_counts,r)

  for (locus = 0; locus < data.nbr_loci; locus++){
    //printf ("Locus = %d of %d\n ",locus+1,data.nbr_loci);

    //declare random number generator and set type to "Mersenne Twister" (MT19937)
    r = gsl_rng_alloc (gsl_rng_mt19937);

    // set seed for random number generator
    //printf ("Random number generator: ’%s’ \n", gsl_rng_name (r));
    gsl_rng_set(r,locus+seed); //change to time so it is different in each run

    one_locus_genotype_counts[0][0] = data.genotype_counts[locus][0][0];
    one_locus_genotype_counts[0][1] = data.genotype_counts[locus][0][1];
    one_locus_genotype_counts[0][2] = data.genotype_counts[locus][0][2];
    one_locus_genotype_counts[1][0] = data.genotype_counts[locus][1][0];
    one_locus_genotype_counts[1][1] = data.genotype_counts[locus][1][1];
    one_locus_genotype_counts[1][2] = data.genotype_counts[locus][1][2];
    if (data.maf[locus] == 1){
      global_result.FST[locus] = one_locus_FST (one_locus_genotype_counts);
      if (user_supplied_fis == 0) {
        if (global_result.Fis < 0){
          global_result.pvalue[locus] = p_value(r, global_result.Ne, 0.0, NBR_SIMULS, one_locus_genotype_counts);
        }else{
          global_result.pvalue[locus] = p_value(r, global_result.Ne, global_result.Fis, NBR_SIMULS, one_locus_genotype_counts);
        }
      } else {
        if (fis < 0){
          global_result.pvalue[locus] = p_value(r, global_result.Ne, 0.0, NBR_SIMULS, one_locus_genotype_counts);
        }else{
          global_result.pvalue[locus] = p_value(r, global_result.Ne, fis, NBR_SIMULS, one_locus_genotype_counts);
        }
      }
      //printf ("Locus %d: Fst = %f; p-value = %f\n", locus+1,global_result.FST[locus],global_result.pvalue[locus]);
    } else {
      global_result.FST[locus] = one_locus_FST (one_locus_genotype_counts);
      global_result.pvalue[locus] = ML_NAN;
      //printf ("Locus %d: Fst = %f; p-value = NA\n", locus+1,global_result.FST[locus]);
      //printf ("Locus %d: Fst = %f; p-value = %f\n", locus+1,global_result.FST[locus],global_result.pvalue[locus]);
    }
    gsl_rng_free (r);

  }

  write_results(&data,&global_result,results_filename);
  write_globalres(&global_result,multilocus_results_filename);


  // need to free memory here?

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
  printf("-maf\t\t\t minumium allele frequeuny threshold for including a locus in the analysis\n");
  printf("    \t\t\t (default 0.0; only values from 0 to 0.5 admited).\n");
  printf("-tau\t\t\t number of generation between the two time samples (no default value)\n");
  printf("-fis\t\t\t inbreeding coefficient, Fis (by default estimated from data)\n");
  printf("-infile\t\t\t input file name (default: \"data.txt\"; must me located in data folder; 100 characters long maximum)\n");
  printf("-outfile\t\t output file name; 100 characters long maximum\n");
  printf("-threads\t\t number of threads, otherwise selected authomatically\n");
  exit(EXIT_SUCCESS);
}

// PRINT VERSION

void print_version()
{
  printf("You are running version %s\n",VERSION);
}
