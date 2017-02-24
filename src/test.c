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

#include "test.h"
#include "readdata.h"
#include "main.h"

void print_test(const gsl_rng * r)
{

  double allele_freq;
  double obs_Fst;
  double p;

  unsigned int one_locus_genotype_counts[2][3];
  unsigned int sample_size[2];
  unsigned int sim_genotype_counts[2][3];

  char *data_filename = "data/data.txt";
  data_struct data;

  // title
  printf("You are running a set of tests to check DriftTest code\n\n");

  // seed
  //printf("Checking RNG seed\n");
  //printf("Seed = %lu\n", seed);
  //if (seed != 123456)
  //{
  //  seed = 123456; 
  //  printf("Seed was reset to %lu\n", seed);
  //}
  //printf("\n");
  
  // drift_sim
  printf("Checking function drift_sim(r,freq,Ne)\n");
  allele_freq = drift_sim(r, 0.5, 10000);
  printf ("Simulated drift with freq=0.5, Ne=10000, tau=10\n");
  printf ("Final allele frequency: %f\n", allele_freq);
  printf ("It should be:  0.4852\n");
  if (allele_freq!= 0.4852) fprintf(stderr, "Warning: final frequency should be 0.4852 for seed=123456\n");
  allele_freq = drift_sim(r, 0.5, 10);
  printf ("Simulated drift with freq=0.5, Ne=10, tau=10\n");
  printf ("Final allele frequency: %f\n", allele_freq);
  printf ("It should be:  0.2\n\n");
  if (allele_freq!= 0.2) fprintf(stderr, "Warning: final frequency should be 0.2 for seed=123456\n\n");

  // counts_sim
  printf("Checking function counts_sim(r,Ne,Fis,one_locus_genotype_counts[],sample_size[],sim_genotype_counts[2][3])\n");
  one_locus_genotype_counts[0][0] = 20;
  one_locus_genotype_counts[0][1] =  5;
  one_locus_genotype_counts[0][2] = 75;
  one_locus_genotype_counts[1][0] = 80;
  one_locus_genotype_counts[1][1] =  2;
  one_locus_genotype_counts[1][2] = 18;
  sample_size[0] = one_locus_genotype_counts[0][0] + one_locus_genotype_counts[0][1] + one_locus_genotype_counts[0][2];
  sample_size[1] = one_locus_genotype_counts[1][0] + one_locus_genotype_counts[1][1] + one_locus_genotype_counts[1][2];
  counts_sim(r, 100, 0.8, one_locus_genotype_counts[0], sample_size, sim_genotype_counts);
  printf("Simulated genotypes with Ne=100, Fis=0.8, initial observed counts 20/5/75:\n");
  printf ("start AA: %d \n", sim_genotype_counts[0][0]);
  printf ("start AB: %d \n", sim_genotype_counts[0][1]);
  printf ("start BB: %d \n", sim_genotype_counts[0][2]);
  printf ("end   AA: %d \n", sim_genotype_counts[1][0]);
  printf ("end   AB: %d \n", sim_genotype_counts[1][1]);
  printf ("end   BB: %d \n", sim_genotype_counts[1][2]);
  printf ("It should be: 16/12/72 & 6/2/92 \n\n");
  if (sim_genotype_counts[0][0]!= 16 || sim_genotype_counts[0][1]!= 12 || sim_genotype_counts[0][2]!= 72)
  {
    fprintf(stderr, "Warning: initial genotype counts should be 16/12/72 for seed=123456\n");
  }
  if (sim_genotype_counts[1][0]!= 6 || sim_genotype_counts[1][1]!= 2 || sim_genotype_counts[1][2]!= 92)
  {
    fprintf(stderr, "Warning: initial genotype counts should be 6/2/92 for seed=123456\n\n");
  }

  // FST_2pop_from_genotypes_counts
  printf("Checking function FST_2pop_from_genotypes_counts(one_locus_genotype_counts[2][3])\n");
  printf("Fst from genotype counts 20/5/75 & 80/2/18:\n");
  obs_Fst = FST_2pop_from_genotypes_counts(one_locus_genotype_counts);
  printf ("Fst = %f\n", obs_Fst);
  printf ("It should be: 0.505721\n");
  //if (obs_Fst!= 0.505721) fprintf(stderr, "Warning: observed Fst should be 0.505721\n");
  printf("Fst from simulated genotype counts (they should be 16/12/72 & 6/2/92):\n");
  obs_Fst = FST_2pop_from_genotypes_counts(sim_genotype_counts);
  printf ("Fst = %f\n", obs_Fst);
  printf ("It should be: 0.078945\n\n");
  //if (obs_Fst!= 0.078945) fprintf(stderr, "Warning: observed Fst should be 0.078945 for seed=123456\n\n");
  
  // p_value
  printf("Checking function p_value(r,Ne,Fis,nbr_simuls,one_locus_genotype_counts[2][3])\n");
  //printf("p-value for genotype counts 20/5/75 & 80/2/18\n  with Ne=100, Fis=0.8 and tau=10 (minimum 1000 simulations):\n");
  //p = p_value(r, 100, 0.8, 1000, one_locus_genotype_counts);
  //printf ("p-value = %f\n", p);
  printf("p-value for simulated genotype counts (they should be 16/12/72 & 6/2/92)\n  with Ne=100, Fis=0.8 and tau=10 (minimum 1000 simulations):\n");
  p = p_value(r, 100, 0.8, 1000, sim_genotype_counts);
  printf ("p-value = %f\n\n", p);

  printf ("Reading data file\n");
  read_data(&data,data_filename);
  
  printf ("Number of individuals in data file = %d\n", data.nbr_ind);
  printf ("Number of loci in data file        = %d\n", data.nbr_loci);
  printf ("First individual belong to population = %d\n", data.pop[0]);
  printf ("data.genotypes[0][0] = %d\n", data.genotypes[0][0]);
  printf ("data.genotypes[0][1] = %d\n", data.genotypes[0][1]);
  printf ("data.genotypes[0][2] = %d\n", data.genotypes[0][2]);
  printf ("data.genotypes[1][0] = %d\n", data.genotypes[1][0]);
  printf ("data.genotypes[1][1] = %d\n", data.genotypes[1][1]);
  printf ("data.genotypes[1][2] = %d\n", data.genotypes[1][2]);


}
