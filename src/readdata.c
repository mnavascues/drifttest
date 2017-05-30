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
  int pop;
  int count,dummy,nloci,nind;
  int min1,max1,min2,max2;
  double mean1,mean2;
  double freq1,freq2,freqMean;

  if ((infile = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"%s: file not found\n",filename);
    exit(EXIT_FAILURE);
  }

  // read first line (number of individuals)
  data -> nbr_ind = 0;
  if ( fscanf(infile,"%d",&data -> nbr_ind) == 1){
    printf("\n\nReading number of individuals: %d \n", data -> nbr_ind);
  } else {
    fprintf(stderr,"Failed to read number of individuals from file (line one)\n");
    check_infile();// program stops
  }
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r')); // read until finding escape character

  // read second line (number of loci)
  data -> nbr_loci = 0;
  if ( fscanf(infile,"%d",&data -> nbr_loci) == 1){
    printf("Reading number of loci: %d \n", data -> nbr_loci);
  } else {
    fprintf(stderr,"Failed to read number of loci from file (line two)\n");
    check_infile();// program stops
  }
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r')); // read until finding escape character

  // allocation memory for genotype counts and sample size
  data -> genotype_counts = (unsigned int ***) malloc(data -> nbr_loci * sizeof(unsigned int**));
  for (i = 0; i < data -> nbr_loci; ++i) {
    data -> genotype_counts[i]=(unsigned int **) malloc(2 * sizeof(unsigned int*));
    for (j = 0; j < 2; ++j) {
      data -> genotype_counts[i][j]=(unsigned int *) malloc(3 * sizeof(unsigned int));
      for (k = 0; k < 3; ++k){
        data -> genotype_counts[i][j][k]=0;
      }
    }
  }
  data -> sample_size = (int **) malloc(data -> nbr_loci * sizeof(int*));
  for (i = 0; i < data -> nbr_loci; ++i) {
    data -> sample_size[i]=(int *) malloc(2 * sizeof(int));
    for (j = 0; j < 2; ++j) {
      data -> sample_size[i][j]=0;
    }
  }
  data -> maf = (int *) malloc(data -> nbr_loci * sizeof(int));


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
          check_infile();// program stops
        }
        if (count != -9 && count != 0 && count != 1 && count != 2) {
          fprintf(stderr,"%s: unexpected count of alleles or time-sample ID at individual %d, column %d.\n",filename,(nind + 1), (nloci + 1));
          fprintf(stderr,"Valid values: -9 for missing data; 0, 1 or 2 copies of reference allele; 1 or 2 for time-sample ID\n");
          check_infile();// program stops
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
            check_infile();// program stops
          }
          if (nloci < (1 + data -> nbr_loci)) {
            fprintf(stderr,"%s: the number of loci at individual no. %d is lower than expected.\n",filename,nind);
            check_infile();// program stops
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
    check_infile();// program stops
  }
  if (nind > data -> nbr_ind) {
    fprintf(stderr,"%s: the number of ind (%d) is larger than expected (%d).\n",filename,nind,data -> nbr_ind);
    check_infile();// program stops
  }
  // printf("OK\n");

  // rewind file and re-read first two lines
  rewind(infile);
  if ( fscanf(infile,"%d",&data -> nbr_ind) != 1){
    fprintf(stderr,"Failed to read number of individuals from file (line one); failed at second pass\n");
    check_infile();// program stops
  }
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
  if ( fscanf(infile,"%d",&data -> nbr_loci) != 1){
     fprintf(stderr,"Failed to read number of loci from file (line two); failed at second pass\n");
     check_infile();// program stops
  }
  while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));

  // read genotypes
  for (j = 0; j < data -> nbr_ind; ++j) {
    // read first value of line : sample of origin of the individual
    if ( fscanf(infile,"%d",&pop) == 1){
      for (i = 0; i < data -> nbr_loci; ++i) {
        if (fscanf(infile,"%d",&count)== 1){
          if (pop == 1){
            if (count == 0) ++data -> genotype_counts[i][0][0];
            if (count == 1) ++data -> genotype_counts[i][0][1];
            if (count == 2) ++data -> genotype_counts[i][0][2];
          }else if (pop == 2){
            if (count == 0) ++data -> genotype_counts[i][1][0];
            if (count == 1) ++data -> genotype_counts[i][1][1];
            if (count == 2) ++data -> genotype_counts[i][1][2];
          }else{
            fprintf(stderr,"Unexpected time-sample ID at individual %d \n", j);
            fprintf(stderr,"Valid values: 1 or 2 (for the recent and ancient sampling times respectively)\n");
            check_infile();// program stops
          }
        } else{
          fprintf(stderr,"Failed to read genotype of locus %d individual %d \n",i,j);
          check_infile();// program stops
        }
      }
    } else {
      fprintf(stderr,"Failed to read time-sample ID (first column) of individual %d \n",j);
      check_infile();// program stops
    }
  }

  fclose(infile);
  min1 = 2147483647;
  max1 = -2147483648;
  mean1 = 0.0;
  min2 = 2147483647;
  max2 = -2147483648;
  mean2 = 0.0;
  for (i = 0; i < data -> nbr_loci; ++i){
    data -> sample_size[i][0]  = data -> genotype_counts[i][0][0] + data -> genotype_counts[i][0][1] + data -> genotype_counts[i][0][2];
    data -> sample_size[i][1]  = data -> genotype_counts[i][1][0] + data -> genotype_counts[i][1][1] + data -> genotype_counts[i][1][2];
    freq1 = ((double)data -> genotype_counts[i][0][1] + 2.0 * (double)data -> genotype_counts[i][0][2]) / (double)data -> sample_size[i][0] / 2.0;
    freq2 = ((double)data -> genotype_counts[i][1][1] + 2.0 * (double)data -> genotype_counts[i][1][2]) / (double)data -> sample_size[i][1] / 2.0;
    freqMean = (freq1+freq2)/2.0;
    if (freqMean <= maf || freqMean >= (1.0 - maf) ){
      data -> maf[i] = 0;
      //printf("Locus %d will not be used in the analysis: minor allele frequency (%f) is lower than the threshold %f\n",i+1,freqMean,maf);
    }else{
      data -> maf[i] = 1;
    }
    mean1 += data -> sample_size[i][0]; 
    mean2 += data -> sample_size[i][1]; 
    if (data -> sample_size[i][0] > max1) max1 = data -> sample_size[i][0];
    if (data -> sample_size[i][0] < min1) min1 = data -> sample_size[i][0];
    if (data -> sample_size[i][1] > max2) max2 = data -> sample_size[i][1];
    if (data -> sample_size[i][1] < min2) min2 = data -> sample_size[i][1];
    //printf("Genotypes locus %d pop1: %d/%d/%d\n",i,data -> genotype_counts[i][0][0],data -> genotype_counts[i][0][1],data -> genotype_counts[i][0][2]);
    //printf("Sample size locus %d pop1: %d\n",i,data -> sample_size[i][0]);
    //printf("Genotypes locus %d pop2: %d/%d/%d\n",i,data -> genotype_counts[i][1][0],data -> genotype_counts[i][1][1],data -> genotype_counts[i][1][2]);
    //printf("Sample size locus %d pop2: %d\n",i,data -> sample_size[i][1]);
  }
  mean1 /= data -> nbr_loci;
  mean2 /= data -> nbr_loci;


  //printf("The data consist in %d individuals typed at %d loci\n",data -> nbr_ind,data -> nbr_loci);
  printf("Mean (min/max) sample size across loci for sample 1; %f (%d,%d)\n",mean1,min1,max1);
  printf("Mean (min/max) sample size across loci for sample 2; %f (%d,%d)\n\n",mean2,min2,max2);

}


void check_infile(){
    fprintf(stderr,"Please check the input file\n");
    exit(EXIT_FAILURE);
}

