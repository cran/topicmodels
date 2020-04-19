#ifndef LDA_INFERENCE_H
#define LDA_INFERENCE_H

#include <math.h>
#include <float.h>
#include <R.h>
#include "lda.h"
#include "utils.h"
#include "common.h"
// #include <assert.h>

double lda_inference(document*, lda_model*, double*, double**, float, int);

#endif
