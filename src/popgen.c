/*
 Copyright INRA
 author: Miguel Navascués (2016)

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

#include "popgen.h"

// simulations of drift using binomial
// ouput: final frequency
double drift_sim(const gsl_rng * r,
                 double freq,
                 unsigned int Ne)
{
  int i,counts;
  
  for (i = 0; i < tau; i++) {
    counts = (int) gsl_ran_binomial (r, freq, Ne);
    freq   = (double) counts / (double) Ne;
  }
  return freq;
}


/* Simulation of genotype frequencies in a sample of idividuals
   taken at the begining and at the end of a simulation of drift
   with function drift_sim()
*/
void counts_sim(const gsl_rng * r,
                unsigned int Ne,
                double Fis,
                const unsigned int genotype_counts[],
                const unsigned int sample_size[],
                unsigned int sim_genotype_counts[2][3])
{
  double alpha[3];
  double genotype_frequencies[3];
  double allele_freq; 

  // alpha parameters for dirichlet
  alpha[0] = 1.0 + (double) genotype_counts[0];
  alpha[1] = 1.0 + (double) genotype_counts[1];
  alpha[2] = 1.0 + (double) genotype_counts[2];

  // use dirichlet-multinomial to get frequencies of genotypes in population and sample
  gsl_ran_dirichlet (r, 3, alpha, genotype_frequencies);
  gsl_ran_multinomial (r, 3, sample_size[0], genotype_frequencies, sim_genotype_counts[0]);
  
  // calculate allele frequency and simulate drift
  allele_freq = genotype_frequencies[0] + genotype_frequencies[1] / 2.0;
  allele_freq = drift_sim (r, allele_freq, Ne);

  // calculate allele frequencies at the end of the simulations
  genotype_frequencies[0] = pow (allele_freq, 2.0) + Fis * (1.0 - allele_freq) * allele_freq;
  genotype_frequencies[1] = 2.0 * (1.0 - allele_freq) * allele_freq * (1.0 - Fis) ;
  genotype_frequencies[2] = pow ( (1.0 - allele_freq), 2.0) + Fis * (1 - allele_freq) * allele_freq;

  gsl_ran_multinomial (r, 3, sample_size[1], genotype_frequencies, sim_genotype_counts[1]);
}



/* single locus genotypes is given as two rows (one for sample)
                and three columns: counts of homozygous allele 1,
                counts of heterozygous and counts of homozygous allele 2    
*/
double FST_2pop_from_genotypes_counts (unsigned int genotype_counts[2][3])
{
  int pop;
  int counts[2] = {0, 0};
  int n[2] = {0, 0};
  int nt = 0;
  int n2 = 0;
  double nc = 0.0;
  double p1[2], p2[2];
  double p1bar, p2bar;
  double frq_hmzgtes_p1[2], frq_hmzgtes_p2[2];
  double SSG = 0;
  double SSI = 0;
  double SSP = 0;
  double MSG, MSI, MSP;
  double FST;

  for (pop = 0; pop < 2; pop++) {
    // get allele counts of allele 2
    counts[pop] = genotype_counts[pop][0] * 2 + genotype_counts[pop][1]; 
    // get sample sizes
    n[pop] = genotype_counts[pop][0] + genotype_counts[pop][1] + genotype_counts[pop][2];
  }

  // get total sample size and other sample size terms for Weir and Cockerham
  nt = n[0] + n[1];
  n2 = (int) pow((double) n[0], 2.0) + (int) pow((double) n[1], 2.0);
  nc = (double) nt - (double) n2 / (double) nt;

  for (pop = 0; pop < 2; pop++) {
    p1[pop] = (2.0 * n[pop] - counts[pop]) / (2.0 * n[pop]);
    p2[pop] = counts[pop] / (2.0 * n[pop]);
 
    frq_hmzgtes_p1[pop] = (double) genotype_counts[pop][0] / (double) n[pop];
    frq_hmzgtes_p2[pop] = (double) genotype_counts[pop][2] / (double) n[pop];
  }

  p1bar = (2.0 * n[0] - counts[0] + 2.0 * n[1] - counts[1]) / (2.0 * n[0] + 2.0 * n[1]);
  p2bar = (counts[0] + counts[1]) / (2.0 * n[0] + 2.0 * n[1]);

  for (pop = 0; pop < 2; pop++) {
    SSG += (double) n[pop] * (p1[pop] - frq_hmzgtes_p1[pop]);
    SSG += (double) n[pop] * (p2[pop] - frq_hmzgtes_p2[pop]);

    SSI += (double) n[pop] * (p1[pop] + frq_hmzgtes_p1[pop] - 2.0 * pow (p1[pop], 2.0) );
    SSI += (double) n[pop] * (p2[pop] + frq_hmzgtes_p2[pop] - 2.0 * pow (p2[pop], 2.0) );

    SSP += (double) n[pop] * pow (p1[pop] - p1bar, 2.0);
    SSP += (double) n[pop] * pow (p2[pop] - p2bar, 2.0);
  }
  SSP *= 2.0;
  MSG = SSG / (double) nt;
  MSI = SSI / ((double) nt - 2.0);
  MSP = SSP;

  FST = (MSP - MSI) / (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);
  //FIS = (MSI - MSG) / (MSI + MSG);
  return FST;
}
                

double p_value(const gsl_rng * r,
               unsigned int Ne,
               double Fis,
               int nbr_simuls,
               unsigned int genotype_counts[2][3])
{
  int pop, sim;
  double obs_Fst, sim_Fst;
  double allele_freq[2], mean_allele_freq;
  unsigned int sample_size[2];
  unsigned int sim_genotype_counts[2][3];
  double p_value = 0.0;

  for (pop = 0; pop < 2; pop++) {
    sample_size[pop] = genotype_counts[pop][0] + genotype_counts[pop][1] + genotype_counts[pop][2];
  }
  obs_Fst = FST_2pop_from_genotypes_counts(genotype_counts);

  sim = 0;
  while (sim < nbr_simuls) {
    counts_sim(r, Ne, Fis, genotype_counts[0], sample_size, sim_genotype_counts);
    for (pop = 0; pop < 2; pop++)
    {
      allele_freq[pop] = (sim_genotype_counts[pop][0] * 2 + sim_genotype_counts[pop][1]) / (2.0 * sample_size[pop]);
    }
    mean_allele_freq = (allele_freq[0] + allele_freq[1]) / 2.0;

    if ((mean_allele_freq >= maf) && (mean_allele_freq <= (1.0 - maf)))
    {
      sim_Fst = FST_2pop_from_genotypes_counts(sim_genotype_counts);
      if (sim_Fst >= obs_Fst) p_value += 1.0;
      sim++;
    }
    if ((sim == nbr_simuls) && (p_value <= 30) && (nbr_simuls <= 1e6)) {
      nbr_simuls *= 10;
    }
  }
  p_value /= (double) nbr_simuls;
  return p_value; 
}



