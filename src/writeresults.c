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

#include "writeresults.h"

void write_results(data_struct *data,
                   global_result_struct *global_result,
                   char *filename)
{

  FILE *outfile = NULL;

  if ((outfile = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"%s: file not found\n",filename);
    exit(EXIT_FAILURE);
  }


  int locus;

  fprintf(outfile,"       locus        F_ST        p_value\n");
  for (locus = 0; locus < data -> nbr_loci; ++locus) {
    if (data -> maf[locus] == 1){
      fprintf(outfile,"%12d%12.6f%15.5e\n",(locus + 1),global_result -> FST[locus],global_result -> pvalue[locus]);
    }else{
      fprintf(outfile,"%12d%12.6f             NA\n",(locus + 1),global_result -> FST[locus]);
    }
  }
  fclose(outfile);
}

