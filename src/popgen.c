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
                const unsigned int one_locus_genotype_counts[],
                const unsigned int sample_size[],
                unsigned int sim_genotype_counts[2][3])
{
  double alpha[3];
  double genotype_frequencies[3];
  double allele_freq;

  // alpha parameters for dirichlet
  alpha[0] = 1.0 + (double) one_locus_genotype_counts[0];
  alpha[1] = 1.0 + (double) one_locus_genotype_counts[1];
  alpha[2] = 1.0 + (double) one_locus_genotype_counts[2];

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
double one_locus_FST (unsigned int one_locus_genotype_counts[2][3])
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
  double num, den;

  num = den = 0.0;
  for (pop = 0; pop < 2; pop++) {
    // get allele counts of allele 2
    counts[pop] = one_locus_genotype_counts[pop][0] * 2 + one_locus_genotype_counts[pop][1]; 
    // get sample sizes
    n[pop] = one_locus_genotype_counts[pop][0] + one_locus_genotype_counts[pop][1] + one_locus_genotype_counts[pop][2];
  }

  // get total sample size and other sample size terms for Weir and Cockerham
  nt = n[0] + n[1];
  n2 = (int) pow((double) n[0], 2.0) + (int) pow((double) n[1], 2.0);
  nc = (double) nt - (double) n2 / (double) nt;

  for (pop = 0; pop < 2; pop++) {
    p1[pop] = (2.0 * n[pop] - counts[pop]) / (2.0 * n[pop]);
    p2[pop] = counts[pop] / (2.0 * n[pop]);
 
    frq_hmzgtes_p1[pop] = (double) one_locus_genotype_counts[pop][0] / (double) n[pop];
    frq_hmzgtes_p2[pop] = (double) one_locus_genotype_counts[pop][2] / (double) n[pop];
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

  // FST = (MSP - MSI) / (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);
  // FIS = (MSI - MSG) / (MSI + MSG);

  num = (MSP - MSI);
  den = (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);

  if (den > 0) {
    return ( num / den);
  }
  else {
    return ML_NAN;
  }
}
                

double p_value(const gsl_rng * r,
               unsigned int Ne,
               double Fis,
               int nbr_simuls,
               unsigned int one_locus_genotype_counts[2][3])
{
  int pop, sim;
  double obs_Fst, sim_Fst;
  double allele_freq[2], mean_allele_freq;
  unsigned int sample_size[2];
  unsigned int sim_genotype_counts[2][3];
  double p_value = 0.0;
  
  //printf ("Number of simulations = %d\n",nbr_simuls);

  for (pop = 0; pop < 2; pop++) {
    sample_size[pop] = one_locus_genotype_counts[pop][0] + one_locus_genotype_counts[pop][1] + one_locus_genotype_counts[pop][2];
  }
  obs_Fst = one_locus_FST(one_locus_genotype_counts);

  sim = 0;
  while (sim < nbr_simuls) {
    counts_sim(r, Ne, Fis, one_locus_genotype_counts[0], sample_size, sim_genotype_counts);
    for (pop = 0; pop < 2; pop++)
    {
      allele_freq[pop] = (sim_genotype_counts[pop][0] * 2 + sim_genotype_counts[pop][1]) / (2.0 * sample_size[pop]);
    }
    mean_allele_freq = (allele_freq[0] + allele_freq[1]) / 2.0;

    if ((mean_allele_freq >= maf) && (mean_allele_freq <= (1.0 - maf)))
    {
      sim_Fst = one_locus_FST(sim_genotype_counts);
      if (sim_Fst >= obs_Fst) p_value += 1.0;
      sim++;
    }
  }
  while ((p_value == 0.0) && (nbr_simuls<1e7)) {
    nbr_simuls = nbr_simuls*10;
    //printf ("Number of simulations increased to = %d\n",nbr_simuls);
    //printf ("Number of simulations done         = %d\n",sim);
    while (sim < nbr_simuls) {
      counts_sim(r, Ne, Fis, one_locus_genotype_counts[0], sample_size, sim_genotype_counts);
      for (pop = 0; pop < 2; pop++)
      {
        allele_freq[pop] = (sim_genotype_counts[pop][0] * 2 + sim_genotype_counts[pop][1]) / (2.0 * sample_size[pop]);
      }
      mean_allele_freq = (allele_freq[0] + allele_freq[1]) / 2.0;

      if ((mean_allele_freq >= maf) && (mean_allele_freq <= (1.0 - maf)))
      {
        sim_Fst = one_locus_FST(sim_genotype_counts);
        if (sim_Fst >= obs_Fst) p_value += 1.0;
        sim++;
      }
    }
  }
  /*
  if (p_value == 0.0){
    nbr_simuls = 1e8;
    p_value = 1;
  } 
  */

  p_value /= (double) nbr_simuls;
  return p_value; 
}


/* multilocus genotype counts are given as a 3D array,
   first dimension:  locus
   second dimension: sample (only two samples, ancestral and modern)
   third dimension:  genotype count (only three values, i.e. biallelic data)
                     - counts of homozygous allele 1
                     - counts of heterozygous
                     - counts of homozygous allele 2    
*/
void F_statistics (data_struct *data, global_result_struct *global_result){

  int locus;
  int pop;
  int counts[2] = {0, 0};
  int n[2] = {0, 0};
  int nt = 0;
  int n2 = 0;
  double nc = 0.0;
  double p1[2], p2[2];
  double p1bar, p2bar;
  double frq_hmzgtes_p1[2], frq_hmzgtes_p2[2];
  double SSG;
  double SSI;
  double SSP;
  double MSG, MSI, MSP;
  double numFst, denFst;
  double numFis, denFis;
  
  numFst = denFst = 0.0;
  numFis = denFis = 0.0;
  for (locus = 0; locus < data -> nbr_loci; locus++){

    //printf ("Locus = %d \n ",locus);


    if (data -> maf[locus] == 1){
      counts[0] = counts[1] = 0;
      n[0] = n[1] = 0;
      nt = n2 = 0;
      nc = 0.0;
 
      for (pop = 0; pop < 2; pop++) {
        // get allele counts of allele 2
        counts[pop] = data -> genotype_counts[locus][pop][0] * 2 + data -> genotype_counts[locus][pop][1]; 
        // get sample sizes
        n[pop] = data -> sample_size[locus][pop];
      }

      // get total sample size and other sample size terms for Weir and Cockerham
      nt = n[0] + n[1];
      n2 = (int) pow((double) n[0], 2.0) + (int) pow((double) n[1], 2.0);
      nc = (double) nt - (double) n2 / (double) nt;

      for (pop = 0; pop < 2; pop++) {
        p1[pop] = (2.0 * n[pop] - counts[pop]) / (2.0 * n[pop]);
        p2[pop] = counts[pop] / (2.0 * n[pop]);
 
        frq_hmzgtes_p1[pop] = (double) data -> genotype_counts[locus][pop][0] / (double) n[pop];
        frq_hmzgtes_p2[pop] = (double) data -> genotype_counts[locus][pop][2] / (double) n[pop];
      }

      p1bar = (2.0 * n[0] - counts[0] + 2.0 * n[1] - counts[1]) / (2.0 * n[0] + 2.0 * n[1]);
      p2bar = (counts[0] + counts[1]) / (2.0 * n[0] + 2.0 * n[1]);

      SSG = SSI = SSP = 0;
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

      // FST = (MSP - MSI) / (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);
      // FIS = (MSI - MSG) / (MSI + MSG);
 
      numFst += (MSP - MSI);
      denFst += (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);

 


   
      numFis += (MSI - MSG);
      denFis += (MSI + MSG);




    }
    //printf ("numFst = %f\n",numFst);
    //printf ("denFst = %f\n",numFst);

    //printf ("numFis = %f\n",numFis);
    //printf ("denFis = %f\n",numFis);
  }

  //printf ("numFst = %f\n",numFst);
  //printf ("denFst = %f\n",numFst);
  //printf ("numFis = %f\n",numFis);
  //printf ("denFis = %f\n",numFis);


  global_result -> Fst = numFst/denFst;
  global_result -> Fis = numFis/denFis;
  global_result -> Ne = (int) tau * (1.0 - global_result -> Fst) / 2.0 / global_result -> Fst;

}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//  BOOTSTRAP FOR CONFIDENCE INTERVALS 

void F_statistics_bootstrap (data_struct *data, global_result_struct *global_result){

  int locus, locus2;
  int pop;
  int counts[2] = {0, 0};
  int n[2] = {0, 0};
  int nt = 0;
  int n2 = 0;
  double nc = 0.0;
  double p1[2], p2[2];
  double p1bar, p2bar;
  double frq_hmzgtes_p1[2], frq_hmzgtes_p2[2];
  double SSG;
  double SSI;
  double SSP;
  double MSG, MSI, MSP;
  double numFst, denFst;
  double numFis, denFis;

  double *VnumFst, *VdenFst;
  double *VnumFis, *VdenFis;
  double epsln = 0.00002;
  double t0_Fst, t0_Fis;
  double *tm_Fst, *tm_Fis;
  double *tp_Fst, *tp_Fis;
  double weight;
  double *dt_Fst, *dt_Fis;
  double *ddt_Fst, *ddt_Fis;

  double alpha = 0.95;
  double lower_bound, upper_bound;

  double sum_dt_Fst_square, sum_dt_Fis_square;
  double sum_dt_Fst_cube, sum_dt_Fis_cube;
  double sigma_hat_Fst, sigma_hat_Fis; 
  double a_hat_Fst, a_hat_Fis;	
  double *delta_Fst, *delta_Fis;
  double tpd_Fst, tmd_Fst, tpd_Fis, tmd_Fis;
  double cq_Fst, b_hat_Fst, sum_ddt_Fst;
  double cq_Fis, b_hat_Fis, sum_ddt_Fis;  

  double tmp_Fst, tmp_Fis;
  double z0_hat_Fst, z0_hat_Fis;
  double w_Fst, lambda_Fst, w_Fis, lambda_Fis;
  double tinf_Fst, tinf_Fis, tsup_Fst, tsup_Fis;



  VnumFst = (double *) malloc(data -> nbr_loci * sizeof(double));
  VdenFst = (double *) malloc(data -> nbr_loci * sizeof(double));
  VnumFis = (double *) malloc(data -> nbr_loci * sizeof(double));
  VdenFis = (double *) malloc(data -> nbr_loci * sizeof(double));
  tm_Fst = (double *) malloc(data -> nbr_loci * sizeof(double));
  tm_Fis = (double *) malloc(data -> nbr_loci * sizeof(double));
  tp_Fst = (double *) malloc(data -> nbr_loci * sizeof(double));
  tp_Fis = (double *) malloc(data -> nbr_loci * sizeof(double));
  dt_Fst = (double *) malloc(data -> nbr_loci * sizeof(double));
  dt_Fis = (double *) malloc(data -> nbr_loci * sizeof(double));
  ddt_Fst = (double *) malloc(data -> nbr_loci * sizeof(double));
  ddt_Fis = (double *) malloc(data -> nbr_loci * sizeof(double));
  delta_Fst = (double *) malloc(data -> nbr_loci * sizeof(double));
  delta_Fis = (double *) malloc(data -> nbr_loci * sizeof(double));


  lower_bound = (1.0 - alpha) / 2.0;
  upper_bound = (1.0 + alpha) / 2.0;

  for (locus = 0; locus < data -> nbr_loci; locus++){

    //printf ("Locus = %d \n ",locus);


    if (data -> maf[locus] == 1){
      counts[0] = counts[1] = 0;
      n[0] = n[1] = 0;
      nt = n2 = 0;
      nc = 0.0;
 
      for (pop = 0; pop < 2; pop++) {
        // get allele counts of allele 2
        counts[pop] = data -> genotype_counts[locus][pop][0] * 2 + data -> genotype_counts[locus][pop][1]; 
        // get sample sizes
        n[pop] = data -> sample_size[locus][pop];
      }

      // get total sample size and other sample size terms for Weir and Cockerham
      nt = n[0] + n[1];
      n2 = (int) pow((double) n[0], 2.0) + (int) pow((double) n[1], 2.0);
      nc = (double) nt - (double) n2 / (double) nt;

      for (pop = 0; pop < 2; pop++) {
        p1[pop] = (2.0 * n[pop] - counts[pop]) / (2.0 * n[pop]);
        p2[pop] = counts[pop] / (2.0 * n[pop]);
 
        frq_hmzgtes_p1[pop] = (double) data -> genotype_counts[locus][pop][0] / (double) n[pop];
        frq_hmzgtes_p2[pop] = (double) data -> genotype_counts[locus][pop][2] / (double) n[pop];
      }

      p1bar = (2.0 * n[0] - counts[0] + 2.0 * n[1] - counts[1]) / (2.0 * n[0] + 2.0 * n[1]);
      p2bar = (counts[0] + counts[1]) / (2.0 * n[0] + 2.0 * n[1]);

      SSG = SSI = SSP = 0;
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

      // FST = (MSP - MSI) / (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);
      // FIS = (MSI - MSG) / (MSI + MSG);
 
      VnumFst[locus] = (MSP - MSI);
      VdenFst[locus] = (MSP + ((double) nc - 1.0) * MSI + (double) nc * MSG);

      VnumFis[locus] = (MSI - MSG);
      VdenFis[locus] = (MSI + MSG);
    }
    //printf ("numFst = %f\n",numFst);
    //printf ("denFst = %f\n",numFst);

    //printf ("numFis = %f\n",numFis);
    //printf ("denFis = %f\n",numFis);
  }

  //printf ("numFst = %f\n",numFst);
  //printf ("denFst = %f\n",numFst);
  //printf ("numFis = %f\n",numFis);
  //printf ("denFis = %f\n",numFis);


  numFst = denFst = 0.0;
  numFis = denFis = 0.0;
  for (locus = 0; locus < data -> nbr_loci; locus++){
    numFst += VnumFst[locus];
    denFst += VdenFst[locus];
    numFis += VnumFis[locus];
    denFis += VdenFis[locus];
  }
  t0_Fst = numFst/denFst;
  t0_Fis = numFis/denFis;
  
  for (locus = 0; locus < data -> nbr_loci; locus++){
    numFst = denFst = 0.0;
    numFis = denFis = 0.0;
    for (locus2 = 0; locus2 < data -> nbr_loci; locus2++){
      weight = epsln / data -> nbr_loci + 1.0 / data -> nbr_loci;
      if (locus==locus2){
        weight =  weight - epsln;
      }
      numFst += VnumFst[locus2] * weight;
      denFst += VdenFst[locus2] * weight;
      numFis += VnumFis[locus2] * weight;
      denFis += VdenFis[locus2] * weight;
    }
    tm_Fst[locus] = numFst/denFst;
    tm_Fis[locus] = numFis/denFis;

    numFst = denFst = 0.0;
    numFis = denFis = 0.0;
    for (locus2 = 0; locus2 < data -> nbr_loci; locus2++){
      weight = -epsln / data -> nbr_loci + 1.0 / data -> nbr_loci;
      if (locus==locus2){
        weight =  weight + epsln;
      }
      numFst += VnumFst[locus2] * weight;
      denFst += VdenFst[locus2] * weight;
      numFis += VnumFis[locus2] * weight;
      denFis += VdenFis[locus2] * weight;
    }
    tp_Fst[locus] = numFst/denFst;
    tp_Fis[locus] = numFis/denFis;

    dt_Fst[locus] = (tp_Fst[locus] - tm_Fst[locus]) / (2 * epsln);
    dt_Fis[locus] = (tp_Fis[locus] - tm_Fis[locus]) / (2 * epsln);

    ddt_Fst[locus] =  (tp_Fst[locus] + tm_Fst[locus] - 2 * t0_Fst) / (epsln*epsln);
    ddt_Fis[locus] =  (tp_Fis[locus] + tm_Fis[locus] - 2 * t0_Fis) / (epsln*epsln);

  }

  sum_dt_Fst_square = 0;
  sum_dt_Fis_square = 0;
  sum_dt_Fst_cube = 0;
  sum_dt_Fis_cube = 0;
  sum_ddt_Fst = 0;
  sum_ddt_Fis = 0;
  for (locus = 0; locus < data -> nbr_loci; locus++){
    sum_dt_Fst_square += dt_Fst[locus] * dt_Fst[locus];
    sum_dt_Fis_square += dt_Fis[locus] * dt_Fis[locus];
    sum_dt_Fst_cube   += dt_Fst[locus] * dt_Fst[locus] * dt_Fst[locus];
    sum_dt_Fis_cube   += dt_Fis[locus] * dt_Fis[locus] * dt_Fis[locus];
    sum_ddt_Fst += ddt_Fst[locus];
    sum_ddt_Fis += ddt_Fis[locus];
  }
  sigma_hat_Fst = 1.0 / data -> nbr_loci * sqrt(sum_dt_Fst_square); 
  sigma_hat_Fis = 1.0 / data -> nbr_loci * sqrt(sum_dt_Fis_square); 
  a_hat_Fst = sum_dt_Fst_cube / (6 * pow(data -> nbr_loci,3) * pow(sigma_hat_Fst,3));	
  a_hat_Fis = sum_dt_Fis_cube / (6 * pow(data -> nbr_loci,3) * pow(sigma_hat_Fis,3));

  for (locus = 0; locus < data -> nbr_loci; locus++){
    delta_Fst[locus] = dt_Fst[locus] / (pow(data -> nbr_loci,2) * sigma_hat_Fst);
    delta_Fis[locus] = dt_Fis[locus] / (pow(data -> nbr_loci,2) * sigma_hat_Fis);
  }	

  numFst = denFst = 0.0;
  numFis = denFis = 0.0;
  for (locus = 0; locus < data -> nbr_loci; locus++){
    weight = delta_Fst[locus] * (-epsln) + 1.0 / data -> nbr_loci;
    numFst += VnumFst[locus] * weight;
    denFst += VdenFst[locus] * weight;
    weight = delta_Fis[locus] * (-epsln) + 1.0 / data -> nbr_loci;
    numFis += VnumFis[locus] * weight;
    denFis += VdenFis[locus] * weight;
  }
  tpd_Fst = numFst/denFst;
  tpd_Fis = numFis/denFis;

  numFst = denFst = 0.0;
  numFis = denFis = 0.0;
  for (locus = 0; locus < data -> nbr_loci; locus++){
    weight = delta_Fst[locus] * epsln + 1.0 / data -> nbr_loci;
    numFst += VnumFst[locus] * weight;
    denFst += VdenFst[locus] * weight;
    weight = delta_Fis[locus] * epsln + 1.0 / data -> nbr_loci;
    numFis += VnumFis[locus] * weight;
    denFis += VdenFis[locus] * weight;
  }
  tmd_Fst = numFst/denFst;
  tmd_Fis = numFis/denFis;

  cq_Fst = (tpd_Fst + tmd_Fst - 2.0 * t0_Fst) / (2.0 * sigma_hat_Fst * pow(epsln,2));
  b_hat_Fst = sum_ddt_Fst / (2 * pow(data -> nbr_loci,2));
  cq_Fis = (tpd_Fis + tmd_Fis - 2.0 * t0_Fis) / (2.0 * sigma_hat_Fis * pow(epsln,2));
  b_hat_Fis = sum_ddt_Fis / (2 * pow(data -> nbr_loci,2));

  tmp_Fst = 2 * gsl_cdf_ugaussian_P(a_hat_Fst) * gsl_cdf_ugaussian_P(cq_Fst - b_hat_Fst / sigma_hat_Fst);
  tmp_Fis = 2 * gsl_cdf_ugaussian_P(a_hat_Fis) * gsl_cdf_ugaussian_P(cq_Fis - b_hat_Fis / sigma_hat_Fis);

  if ((tmp_Fst >= 1.0) || (tmp_Fst <= 0.0)) {
    printf("Confidence interval for Fst has infinite range...\n");
  } else {
    z0_hat_Fst = gsl_cdf_ugaussian_Pinv(tmp_Fst);
  }

  if ((tmp_Fis >= 1.0) || (tmp_Fis <= 0.0)) {
    printf("Confidence interval for Fis has infinite range...\n");
  } else {
    z0_hat_Fis = gsl_cdf_ugaussian_Pinv(tmp_Fis);
  }

  w_Fst = z0_hat_Fst + gsl_cdf_ugaussian_Pinv(lower_bound);	  
  lambda_Fst = w_Fst / pow(1.0 - a_hat_Fst * w_Fst, 2);
  w_Fis = z0_hat_Fis + gsl_cdf_ugaussian_Pinv(lower_bound);	  
  lambda_Fis = w_Fis / pow(1.0 - a_hat_Fis * w_Fis, 2);

  numFst = denFst = 0.0;
  numFis = denFis = 0.0;
  for (locus = 0; locus < data -> nbr_loci; locus++){
    weight = delta_Fst[locus] * lambda_Fst + 1.0 / data -> nbr_loci;
    numFst += VnumFst[locus] * weight;
    denFst += VdenFst[locus] * weight;
    weight = delta_Fis[locus] * lambda_Fis + 1.0 / data -> nbr_loci;
    numFis += VnumFis[locus] * weight;
    denFis += VdenFis[locus] * weight;
  }
  tinf_Fst = numFst/denFst;
  tinf_Fis = numFis/denFis;


  w_Fst = z0_hat_Fst + gsl_cdf_ugaussian_Pinv(upper_bound);	  
  lambda_Fst = w_Fst / pow(1.0 - a_hat_Fst * w_Fst, 2);
  w_Fis = z0_hat_Fis + gsl_cdf_ugaussian_Pinv(upper_bound);	  
  lambda_Fis = w_Fis / pow(1.0 - a_hat_Fis * w_Fis, 2);

  numFst = denFst = 0.0;
  numFis = denFis = 0.0;
  for (locus = 0; locus < data -> nbr_loci; locus++){
    weight = delta_Fst[locus] * lambda_Fst + 1.0 / data -> nbr_loci;
    numFst += VnumFst[locus] * weight;
    denFst += VdenFst[locus] * weight;
    weight = delta_Fis[locus] * lambda_Fis + 1.0 / data -> nbr_loci;
    numFis += VnumFis[locus] * weight;
    denFis += VdenFis[locus] * weight;
  }
  tsup_Fst = numFst/denFst;
  tsup_Fis = numFis/denFis;


  global_result -> Fst = t0_Fst;
  global_result -> Fst_lower = tinf_Fst;
  global_result -> Fst_upper = tsup_Fst;


  global_result -> Fis = t0_Fis;
  global_result -> Fis_lower = tinf_Fis;
  global_result -> Fis_upper = tsup_Fis;

  global_result -> Ne = (int) tau * (1.0 - t0_Fst) / 2.0 / t0_Fst;
  global_result -> Ne_lower = (int) tau * (1.0 - tsup_Fst) / 2.0 / tsup_Fst;
  global_result -> Ne_upper = (int) tau * (1.0 - tinf_Fst) / 2.0 / tinf_Fst;

  global_result -> sigma = 2*t0_Fis / (1+t0_Fis);
  global_result -> sigma_lower = 2*tinf_Fis / (1+tinf_Fis);
  global_result -> sigma_upper = 2*tsup_Fis / (1+tsup_Fis);


}

