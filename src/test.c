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

void print_test(const gsl_rng * r)
{

  printf("You are running a set of tests to verify the integrity of DriftTest code\n\n");

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
