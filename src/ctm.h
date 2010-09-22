// (C) Copyright 2007, David M. Blei and John D. Lafferty

// This file is part of CTM-C.

// CTM-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// CTM-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA


#ifndef LLNA_H
#define LLNA_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <time.h>
#include <R.h>

#define NUM_INIT 1
#define SEED_INIT_SMOOTH 1.0

/*
 * the llna model
 *
 */

typedef struct llna_model
{
    int k;
    gsl_matrix * log_beta;
    gsl_vector * mu;
    gsl_matrix * inv_cov;
    gsl_matrix * cov;
    double log_det_inv_cov;
} llna_model;


/*
 * sufficient statistics for mle of an llna model
 *
 */

typedef struct llna_ss
{
    gsl_matrix * cov_ss;
    gsl_vector * mu_ss;
    gsl_matrix * beta_ss;
    double ndata;
} llna_ss;


/*
 * BG: llna_params taken from params.h
 *
 */

typedef struct llna_params
{
    int em_max_iter;
    int var_max_iter;
    int cg_max_iter;
    double em_convergence;
    double var_convergence;
    double cg_convergence;
    int cov_estimate;
    int verbose;
    int save;
    int seed;
} llna_params;

/*
 * BG: doc and corpus taken from corpus.h
 * a document is a collection of counts and terms
 *
 */

typedef struct doc {
    int total;
    int nterms;
    int * word;
    int * count;
} doc;


/*
 * a corpus is a collection of documents
 *
 */

typedef struct corpus {
    doc* docs;
    int nterms;
    int ndocs;
} corpus;



/*
 * function declarations
 *
 */

void write_llna_model(llna_model*, char*, int);
llna_model* new_llna_model(int, int);
llna_model* random_init(int, int, int, int);
llna_model* corpus_init(int, corpus*, int, int);
llna_ss * new_llna_ss(llna_model*);
void del_llna_ss(llna_ss*);
void reset_llna_ss(llna_ss*);
void write_ss(llna_ss*);

#endif
