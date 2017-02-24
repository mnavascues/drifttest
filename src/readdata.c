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

#include "readdata.h"

void read_data(data_struct *data,
               char *filename)
{
  FILE *infile = NULL;
  char X;
  int i,j,k;
  int count,dummy,nloci,nind;

  if ((infile = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"%s: file not found\n",filename);
    exit(EXIT_FAILURE);
  }

  // read first line (number of individuals)
  data -> nbr_ind = 0;
  fscanf(infile,"%d",&data -> nbr_ind);
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r')); // read until finding escape character

  // read second line (number of loci)
  data -> nbr_loci = 0;
  fscanf(infile,"%d",&data -> nbr_loci);
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r')); // read until finding escape character

  // allocation momory for genotypes and population origin
  data -> genotypes = (int **) malloc(data -> nbr_loci * sizeof(int *));
  for (i = 0; i < data -> nbr_loci; ++i) {
    data -> genotypes[i]=(int *) malloc(data -> nbr_ind * sizeof(int));
    for (j = 0; j < data -> nbr_ind; ++j) {
      data -> genotypes[i][j] = -99;
    }
  }
  data -> pop = (int *) malloc(data -> nbr_ind * sizeof(int));
  for (j = 0; j < data -> nbr_ind; ++j) {
    data -> pop[j] = 0;
  }
  data -> genotype_counts = (int ***) malloc(data -> nbr_loci * sizeof(int**));
  for (i = 0; i < data -> nbr_loci; ++i) {
    data -> genotype_counts[i]=(int **) malloc(2 * sizeof(int*));
    for (j = 0; j < 2; ++j) {
      data -> genotype_counts[i][j]=(int *) malloc(3 * sizeof(int));
      for (k = 0; k < 3; ++k){
        data -> genotype_counts[i][j][k]=0;
      }
    }
  }

  //verification number of loci and individuals
  nind = 0;
  do {
    nloci = 0;
    do {
      dummy = -1;
      X = fscanf(infile,"%i%n",&count,&dummy);
      //printf("X=%d; count=%d; dummy=%d\n",X,count,dummy);
      if (X != EOF) {
        if (dummy == -1) {
          fprintf(stderr,"%s: unexpected character at individual %d, locus %d.\n",filename,(nind + 1), (nloci + 1));
          check_infile();
        }
        if (count != -9 && count != 0 && count != 1 && count != 2) {
          fprintf(stderr,"%s: unexpected count of alleles at individual %d, locus %d.\n",filename,(nind + 1), (nloci + 1));
          fprintf(stderr,"Valid values: -9, 0, 1, 2.\n");
          check_infile();
        }
      }
      nloci++;
      do {
        X = getc(infile);
	//printf("X=%d\n",X);
        if (X == '\n') {
          nind++;
          if (nloci > (1 + data -> nbr_loci)) {
            fprintf(stderr,"%s: the number of loci at individual no. %d is larger than expected.\n",filename,nind);
            check_infile();
          }
          if (nloci < (1 + data -> nbr_loci)) {
            fprintf(stderr,"%s: the number of loci at individual no. %d is lower than expected.\n",filename,nind);
            check_infile();
          }
          if (nloci == (1 + data -> nbr_loci)) {
            break;
          }
        }
      } while (X  == '\t' || X == ' ');
      fseek(infile, -1, SEEK_CUR);
    } while ((X != '\n') & (X != EOF));
  } while (X != EOF);
  if (nind < data -> nbr_ind) {
    fprintf(stderr,"%s: the number of ind (%d) is lower than expected (%d).\n",filename,nind,data -> nbr_ind);
    check_infile();
  }
  if (nind > data -> nbr_ind) {
    fprintf(stderr,"%s: the number of ind (%d) is larger than expected (%d).\n",filename,nind,data -> nbr_ind);
    check_infile();
  }
  printf("OK\n");

  // rewind file and re-read first two lines
  rewind(infile);
  fscanf(infile,"%d",&data -> nbr_ind);
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
  fscanf(infile,"%d",&data -> nbr_loci);
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));

  // read genotypes
  for (j = 0; j < data -> nbr_ind; ++j) {
    fscanf(infile,"%d",&count);
    data -> pop[j] = count; 
    for (i = 0; i < data -> nbr_loci; ++i) {
      fscanf(infile,"%d",&count);
      data -> genotypes[j][i] = count; 
    }
  }

  fclose(infile);
  printf("The data consist in %d individuals typed at %d loci\n\n",data -> nbr_ind,data -> nbr_loci);
  
  for (i = 0; i < data -> nbr_loci; ++i){
    for (j = 0; j < data -> nbr_ind; ++j){
      if (data -> pop[j] == 1){
        if (data -> genotypes[j][i] == 0) ++data -> genotype_counts[i][0][0];
        if (data -> genotypes[j][i] == 1) ++data -> genotype_counts[i][0][1];
        if (data -> genotypes[j][i] == 2) ++data -> genotype_counts[i][0][2];
      }else if (data -> pop[j] == 2){
        if (data -> genotypes[j][i] == 0) ++data -> genotype_counts[i][1][0];
        if (data -> genotypes[j][i] == 1) ++data -> genotype_counts[i][1][1];
        if (data -> genotypes[j][i] == 2) ++data -> genotype_counts[i][1][2];
      }
    }
    printf("Genotypes locus %d pop1: %d/%d/%d\n",i,data -> genotype_counts[i][0][0],data -> genotype_counts[i][0][1],data -> genotype_counts[i][0][2]);
    printf("Genotypes locus %d pop2: %d/%d/%d\n",i,data -> genotype_counts[i][1][0],data -> genotype_counts[i][1][1],data -> genotype_counts[i][1][2]);
  }
  


}


void check_infile(){
    fprintf(stderr,"Please check the input file\n");
    exit(EXIT_FAILURE);
}

