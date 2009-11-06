#ifndef RLDA_H
#define RLDA_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>

#include <R.h>
#include <Rdefines.h>
#include "lda.h"
#include "lda-inference.h"
#include "lda-model.h"
#include "lda-alpha.h"
#include "utils.h"
#include "cokus.h"

float EM_CONVERGED;
int EM_MAX_ITER;
int ESTIMATE_ALPHA;
double INITIAL_ALPHA;
int LAG;
int NTOPICS;

SEXP rlda(SEXP i, SEXP j, SEXP v, SEXP nrow, SEXP ncol,
	  SEXP control, SEXP k, SEXP prefix, SEXP init_model);

SEXP rlda_inf(SEXP x, SEXP init_model, SEXP prefix);

#endif
