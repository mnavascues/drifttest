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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>            // for computing common mathematical operations
#include <gsl/gsl_rng.h>     // GSL: random number generation
#include <gsl/gsl_randist.h> // GSL: random number distributions

#include "drifttest.h"

double drift_sim(const gsl_rng * r, double freq, unsigned int Ne);

void counts_sim(const gsl_rng * r, unsigned int Ne,
                double Fis, const unsigned int genotype_counts[],
                const unsigned int sample_size[],
                unsigned int sim_genotype_counts[2][3]);

double FST_2pop_from_genotypes_counts (unsigned int genotype_counts[2][3]);

double p_value(const gsl_rng * r, unsigned int Ne, double Fis,
               int nbr_simuls, unsigned int genotype_counts[2][3]);
