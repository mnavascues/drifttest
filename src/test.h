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
#ifndef _TEST_H
#define _TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>            // for computing common mathematical operations
#include <gsl/gsl_rng.h>     // GSL: random number generation
#include <gsl/gsl_randist.h> // GSL: random number distributions

#include "main.h"
#include "popgen.h"

void print_test(const gsl_rng * r);

#endif
