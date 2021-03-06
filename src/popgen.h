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
#ifndef _POPGEN_H
#define _POPGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>            // for computing common mathematical operations
#include <gsl/gsl_rng.h>     // GSL: random number generation
#include <gsl/gsl_randist.h> // GSL: random number distributions
#include <gsl/gsl_cdf.h>     // GSL: for cdf for probability distributions

#include "main.h"

#define ML_NAN   (0.0 / 0.0) 

double drift_sim(const gsl_rng * r, double freq, unsigned int Ne);

void counts_sim(const gsl_rng * r, unsigned int Ne,
                double Fis, const unsigned int one_locus_genotype_counts[],
                const unsigned int sample_size[],
                unsigned int sim_genotype_counts[2][3]);

double one_locus_FST (unsigned int one_locus_genotype_counts[2][3]);

double p_value(const gsl_rng * r, unsigned int Ne, double Fis,
               int nbr_simuls, unsigned int one_locus_genotype_counts[2][3]);

void F_statistics (data_struct *data, global_result_struct *global_result);
void F_statistics_bootstrap (data_struct *data, global_result_struct *global_result);
#endif
