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
#ifndef _MAIN_H
#define _MAIN_H

#define GLOBAL_VARIABLE

#define VERSION     "1.0.4"
//#define LOGFILE     "results.log"
//#define OUTFILE     "results.out"
#define NBR_SIMULS  1e4
#define ML_NAN      (0.0 / 0.0)

// for text color
#define RED   "\x1B[31m"
#define RESET "\x1B[0m"

// ---------------------
// structure definitions
// ---------------------


typedef struct {
  int nbr_ind;                     // Number of individuals
  int nbr_loci;                    // Number of loci
  int **sample_size;               // sample size per locus per time-sample (after missing data)
  unsigned int ***genotype_counts; // counts of each genotype per locus per time-sample
  int *maf;                        // 0:locus maf<threshold; 1:locus maf>=threshold
} data_struct;

typedef struct {
  double Fst;
  double Fis;
  int Ne;
  double *FST;
  double *pvalue;
} global_result_struct;


// -----------------
//  local functions
// -----------------

void print_usage();
void print_version();

// ------------------
//  GLOBAL VARIABLES
// ------------------

GLOBAL_VARIABLE const char *program_name;
GLOBAL_VARIABLE int tau;
GLOBAL_VARIABLE double maf;
GLOBAL_VARIABLE double fis;

#endif
