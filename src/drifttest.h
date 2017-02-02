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

#define GLOBAL_VARIABLE

#define VERSION     "0.0.0"
#define LOGFILE     "results.log"
#define OUTFILE     "results.out"
#define NBR_SIMULS  1e4

// ---------------------
// structure definitions
// ---------------------

typedef struct {
	int nbr_loci;     // Number of loci
        int *pop;         // population origin of each individual
	int **genotypes;  // genotype for each individual and locus, 
                          // coded as the number of copies of reference allele
} data_struct;


// ----------------
// local functions 
// ----------------

void print_usage();
void print_version();

// ----------------
// GLOBAL VARIABLES
// ----------------

GLOBAL_VARIABLE const char *program_name;
GLOBAL_VARIABLE int tau;
GLOBAL_VARIABLE double maf;
GLOBAL_VARIABLE FILE *logfile;

