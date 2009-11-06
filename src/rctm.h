#ifndef RLDA_H
#define RLDA_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#include <R.h>
#include <Rdefines.h>
#include "ctm.h"
#include "inference.h"
#include "gsl-wrappers.h"

llna_params PARAMS;

SEXP rctm(SEXP i, SEXP j, SEXP v, SEXP nrow, SEXP ncol,
	  SEXP control, SEXP k, SEXP prefix, SEXP init_model);

SEXP rctm_inf(SEXP x, SEXP model_init, SEXP prefix);

#endif
