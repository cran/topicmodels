// (C) Copyright 2007, David M. Blei and John D. Lafferty

// This file is part of CTM-C.

// CTM-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// CTM-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PATICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

/*************************************************************************
 *
 * llna.c
 *
 * reading, writing, and initializing a logistic normal allocation model
 *
 *************************************************************************/

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "gsl-wrappers.h"
#include "ctm.h"

/*
 * create a new empty model and delete a model
 *
 */

llna_model* new_llna_model(int ntopics, int nterms)
{
    llna_model* model = malloc(sizeof(llna_model));
    model->k = ntopics;
    model->mu = gsl_vector_calloc(ntopics - 1);
    model->cov = gsl_matrix_calloc(ntopics-1, ntopics-1);
    model->inv_cov = gsl_matrix_calloc(ntopics-1, ntopics-1);
    model->log_beta = gsl_matrix_calloc(ntopics, nterms);
    return(model);
}

void del_llna_model(llna_model * model)
{
    gsl_vector_free(model->mu);
    gsl_matrix_free(model->cov);
    gsl_matrix_free(model->inv_cov);
    gsl_matrix_free(model->log_beta);
}

/*
 * create and delete sufficient statistics
 *
 */

llna_ss * new_llna_ss(llna_model* model)
{
    llna_ss * ss;
    ss = malloc(sizeof(llna_ss));
    ss->mu_ss = gsl_vector_calloc(model->k-1);
    ss->cov_ss = gsl_matrix_calloc(model->k-1, model->k-1);
    ss->beta_ss = gsl_matrix_calloc(model->k, model->log_beta->size2);
    ss->ndata = 0;
    reset_llna_ss(ss);
    return(ss);
}


void del_llna_ss(llna_ss * ss)
{
    gsl_vector_free(ss->mu_ss);
    gsl_matrix_free(ss->cov_ss);
    gsl_matrix_free(ss->beta_ss);
}

void reset_llna_ss(llna_ss * ss)
{
    gsl_matrix_set_all(ss->beta_ss, 0);
    gsl_matrix_set_all(ss->cov_ss, 0);
    gsl_vector_set_all(ss->mu_ss, 0);
    ss->ndata = 0;
}


/*
 * initialize a model with zero-mean, diagonal covariance gaussian and
 * topics seeded from the corpus
 *
 */

llna_model* corpus_init(int ntopics, corpus* corpus, int verbose, int seed)
{
    llna_model* model = new_llna_model(ntopics, corpus->nterms);
    gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
    doc* doc;
    int i, k, n, d;
    double sum;
    if (verbose > 0) Rprintf("USING %d\n", seed);
    gsl_rng_set(r, seed);

    // gaussian
    for (i = 0; i < ntopics-1; i++)
    {
        vset(model->mu, i, 0);
        mset(model->cov, i, i, 1.0);
    }
    matrix_inverse(model->cov, model->inv_cov);
    model->log_det_inv_cov = log_det(model->inv_cov);

    // topics
    for (k = 0; k < ntopics; k++)
    {
        sum = 0;
        // seed
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(gsl_rng_uniform(r)*corpus->ndocs);
            if (verbose > 0) Rprintf("initialized with document %d\n", d);
            doc = &(corpus->docs[d]);
            for (n = 0; n < doc->nterms; n++)
            {
                minc(model->log_beta, k, doc->word[n], (double) doc->count[n]);
            }
        }
        // smooth
        for (n = 0; n < model->log_beta->size2; n++)
        {
            minc(model->log_beta, k, n, SEED_INIT_SMOOTH + gsl_rng_uniform(r));
            // minc(model->log_beta, k, n, SEED_INIT_SMOOTH);
            sum += mget(model->log_beta, k, n);
        }
        sum = safe_log(sum);
        // normalize
        for (n = 0; n < model->log_beta->size2; n++)
        {
            mset(model->log_beta, k, n,
                 safe_log(mget(model->log_beta, k, n)) - sum);
        }
    }
    gsl_rng_free(r);
    return(model);
}

/*
 * random initialization means zero-mean, diagonal covariance gaussian
 * and randomly generated topics
 *
 */

llna_model* random_init(int ntopics, int nterms, int verbose, int seed)
{
    int i, j;
    double sum, val;
    llna_model* model = new_llna_model(ntopics, nterms);
    gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
    if (verbose > 0) Rprintf("RANDOM SEED = %ld\n", seed);
    gsl_rng_set(r, seed);

    for (i = 0; i < ntopics-1; i++)
    {
        vset(model->mu, i, 0);
        mset(model->cov, i, i, 1.0);
    }
    for (i = 0; i < ntopics; i++)
    {
        sum = 0;
        for (j = 0; j < nterms; j++)
        {
            val = gsl_rng_uniform(r) + 1.0/100;
            sum += val;
            mset(model->log_beta, i, j, val);
        }
        for (j = 0; j < nterms; j++)
            mset(model->log_beta, i, j, log(mget(model->log_beta, i, j) / sum));
    }
    matrix_inverse(model->cov, model->inv_cov);
    model->log_det_inv_cov = log_det(model->inv_cov);

    gsl_rng_free(r);
    return(model);
}

/*
 * write a model
 *
 */

void write_llna_model(llna_model * model, char * root, int verbose)
{
    char filename[260];
    FILE* fileptr;

    // write parameters
    if (verbose > 0) Rprintf("writing params\n");
    snprintf(filename, 260, "%s-param.txt", root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->k);
    fprintf(fileptr, "num_terms %d\n", (int) model->log_beta->size2);
    fclose(fileptr);
    // write gaussian
    if (verbose > 0) Rprintf("writing gaussian\n");
    snprintf(filename, 260, "%s-mu.dat", root);
    printf_vector(filename, model->mu);
    snprintf(filename, 260, "%s-cov.dat", root);
    printf_matrix(filename, model->cov);
    snprintf(filename, 260, "%s-inv-cov.dat", root);
    printf_matrix(filename, model->inv_cov);
    snprintf(filename, 260, "%s-log-det-inv-cov.dat", root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "%lf\n", model->log_det_inv_cov);
    fclose(fileptr);
    // write topic matrix
    if (verbose > 0) Rprintf("writing topics\n");
    snprintf(filename, 260, "%s-log-beta.dat", root);
    printf_matrix(filename, model->log_beta);
}
