#ifndef LDA_ALPHA_H
#define LDA_ALPHA_H

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include "lda.h"
#include "utils.h"

#define NEWTON_THRESH 1e-5
#define MAX_ALPHA_ITER 1000

double opt_alpha(double ss, int D, int K);

#endif
