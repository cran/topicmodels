// (C) Copyright 2007, David M. Blei and John D. Lafferty
// (C) Copyright 2009, Bettina Gruen (Bettina [dot] Gruen [at] wu-wien [dot] ac [dot] at)

// This file uses part of CTM-C.

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

#include "rctm.h"

/*
 * BG: construct return object
 *
 */

SEXP returnObjectCTM(SEXP ans, llna_model* model, corpus* corpus, gsl_matrix* corpus_lambda, gsl_matrix* corpus_nu,
		     gsl_matrix* corpus_phi_sum, llna_var_param** var, gsl_vector* likelihood, int iter, double* logLiks, int keep_iter) 
{
  SEXP tp, I, J, V, wordassign, nms, dn;
  SEXP ROWNAMES, COLNAMES;
  doc* doc;
  int i, j, d, total, n;
  double *m;
  gsl_vector phi_row;

  tp = PROTECT(allocVector(INTSXP, 1));
  *INTEGER(tp) = iter;
  SET_SLOT(ans, install("iter"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocVector(INTSXP, 1));
  *INTEGER(tp) = model->k;
  SET_SLOT(ans, install("k"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocVector(INTSXP, 2));
  INTEGER(tp)[0] = corpus->ndocs;
  INTEGER(tp)[1] = corpus->nterms;
  SET_SLOT(ans, install("Dim"), tp);
  UNPROTECT(1);
  
  tp = PROTECT(allocVector(REALSXP, corpus->ndocs));
  for (i = 0; i < corpus->ndocs; i++) 
    REAL(tp)[i] = gsl_vector_get(likelihood, i);
  SET_SLOT(ans, install("loglikelihood"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocMatrix(REALSXP, model->k, corpus->nterms));
  m  = REAL(tp);
  for (i = 0; i < model->k; i++)
    for (j = 0; j < corpus->nterms; j++)
      m[i + model->k * j] = gsl_matrix_get(model->log_beta, i, j);
  SET_SLOT(ans, install("beta"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocVector(REALSXP, model->k - 1));
  for (i = 0; i < (model->k - 1); i++)
    REAL(tp)[i] = gsl_vector_get(model->mu, i);
  SET_SLOT(ans, install("mu"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocMatrix(REALSXP, model->k - 1, model->k - 1));
  m = REAL(tp);
  for (i = 0; i < (model->k - 1); i++)
    for (j = 0; j < (model->k - 1); j++)
      m[i + (model->k - 1) * j] = gsl_matrix_get(model->cov, i, j);
  SET_SLOT(ans, install("Sigma"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocMatrix(REALSXP, corpus->ndocs, model->k - 1));
  m = REAL(tp);
  for (i = 0; i < corpus->ndocs; i++)
    for (j = 0; j < (model->k-1); j++)
      m[i + corpus->ndocs * j] = gsl_matrix_get(corpus_nu, i, j);
  SET_SLOT(ans, install("nusquared"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocMatrix(REALSXP, corpus->ndocs, model->k - 1));
  m = REAL(tp);
  for (i = 0; i < corpus->ndocs; i++)
    for (j = 0; j < model->k - 1; j++)
      m[i + corpus->ndocs * j] = gsl_matrix_get(corpus_lambda, i, j);
  SET_SLOT(ans, install("gamma"), tp);
  UNPROTECT(1);

  if ((PARAMS.keep > 0) && (PARAMS.em_max_iter > 0)) {
    tp = PROTECT(allocVector(REALSXP, keep_iter));
    for (i = 0; i < keep_iter; i++) 
      REAL(tp)[i] = logLiks[i];
    SET_SLOT(ans, install("logLiks"), tp);
    UNPROTECT(1);
  }

  wordassign = PROTECT(allocVector(VECSXP, 6));
  total = 0;
  for (d = 0; d < corpus->ndocs; d++) {
    doc = &(corpus->docs[d]);
    total += doc->nterms;
  }
  
  I = PROTECT(allocVector(INTSXP, total));
  J = PROTECT(allocVector(INTSXP, total));
  V = PROTECT(allocVector(REALSXP, total));

  INTEGER(I)[0] = 0;
  i = 0;
  for (d = 0; d < corpus->ndocs; d++) {
    doc = &(corpus->docs[d]);
    for (n = 0; n < doc->nterms; n++) {
      INTEGER(I)[i] = d+1;
      INTEGER(J)[i] = doc->word[n] + 1;
      phi_row = gsl_matrix_row(var[d]->phi, n).vector;
      REAL(V)[i] = argmax(&phi_row) + 1;
      i++;
    }
  }

  SET_VECTOR_ELT(wordassign, 0, I); 
  SET_VECTOR_ELT(wordassign, 1, J);
  SET_VECTOR_ELT(wordassign, 2, V);
  UNPROTECT(3);

  tp = PROTECT(allocVector(INTSXP, 1));
  INTEGER(tp)[0] = corpus->ndocs;
  SET_VECTOR_ELT(wordassign, 3, tp);
  UNPROTECT(1);
  tp = PROTECT(allocVector(INTSXP, 1));
  INTEGER(tp)[0] = corpus->nterms;
  SET_VECTOR_ELT(wordassign, 4, tp);
  UNPROTECT(1);

  dn = PROTECT(allocVector(VECSXP, 2));
  ROWNAMES = PROTECT(allocVector(INTSXP, corpus->ndocs));
  COLNAMES = PROTECT(allocVector(INTSXP, corpus->nterms));
  for (d = 0; d < corpus->ndocs; d++) INTEGER(ROWNAMES)[d] = d+1;
  for (d = 0; d < corpus->nterms; d++) INTEGER(COLNAMES)[d] = d+1;
  SET_VECTOR_ELT(dn, 0, AS_CHARACTER(ROWNAMES));
  SET_VECTOR_ELT(dn, 1, AS_CHARACTER(COLNAMES));
  setAttrib(dn, R_NamesSymbol, nms = allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, mkChar("Docs"));
  SET_STRING_ELT(nms, 1, mkChar("Terms"));
  SET_VECTOR_ELT(wordassign, 5, dn);
  UNPROTECT(3);
  
  setAttrib(wordassign, R_NamesSymbol, nms = allocVector(STRSXP, 6));
  SET_STRING_ELT(nms, 0, mkChar("i"));
  SET_STRING_ELT(nms, 1, mkChar("j"));
  SET_STRING_ELT(nms, 2, mkChar("v"));
  SET_STRING_ELT(nms, 3, mkChar("nrow"));
  SET_STRING_ELT(nms, 4, mkChar("ncol"));
  SET_STRING_ELT(nms, 5, mkChar("dimnames"));
  setAttrib(wordassign, R_ClassSymbol, mkString("simple_triplet_matrix"));
  
  SET_SLOT(ans, install("wordassignments"), wordassign);
  UNPROTECT(1);
  return(ans);
}

llna_model* R2llna_model(SEXP init_model) 
{
  llna_model* model;
  int ntopics, nterms, i, j;

  nterms = INTEGER(GET_SLOT(init_model, install("Dim")))[1];
  ntopics = *INTEGER(GET_SLOT(init_model, install("k")));


  // allocate model
  model = new_llna_model(ntopics, nterms);
  
  for (i = 0; i < (ntopics - 1); i++) 
    gsl_vector_set(model->mu, i, REAL(GET_SLOT(init_model, install("mu")))[i]);
  for (i = 0; i < (ntopics - 1); i++)
    for (j = 0; j < (ntopics - 1); j++)
      gsl_matrix_set(model->cov, i, j, REAL(GET_SLOT(init_model, install("Sigma")))[i + (ntopics - 1) * j]);
  matrix_inverse(model->cov, model->inv_cov);
  model->log_det_inv_cov = log_det(model->inv_cov);
  for (i = 0; i < ntopics; i++)
    for (j = 0; j < nterms; j++)
      gsl_matrix_set(model->log_beta, i, j, REAL(GET_SLOT(init_model, install("beta")))[i + ntopics * j]);
  return(model);
}

/*
 * BG: copied from estimate.c
 * e step
 *
 */

void expectation(corpus* corpus, llna_model* model, llna_ss* ss,
                 double* avg_niter, double* total_lhood,
		 gsl_vector* likelihood,
                 gsl_matrix* corpus_lambda, gsl_matrix* corpus_nu,
                 gsl_matrix* corpus_phi_sum,
                 short reset_var, double* converged_pct, llna_var_param** var, int verbose)
{
    int i;
    doc doc;
    double total;
    gsl_vector lambda, nu;
    gsl_vector* phi_sum;

    *avg_niter = 0.0;
    *converged_pct = 0;
    phi_sum = gsl_vector_alloc(model->k);
    total = 0;
    for (i = 0; i < corpus->ndocs; i++)
    {
      if ((verbose > 0) && ((i % (corpus->ndocs-1)) == 0) && (i>0)) Rprintf("doc %5d   ", i);
        doc = corpus->docs[i];
        var[i] = new_llna_var_param(doc.nterms, model->k);
        if (reset_var)
            init_var_unif(var[i], &doc, model);
        else
        {
            lambda = gsl_matrix_row(corpus_lambda, i).vector;
            nu= gsl_matrix_row(corpus_nu, i).vector;
            init_var(var[i], &doc, model, &lambda, &nu);
        }
        gsl_vector_set(likelihood, i, var_inference(var[i], &doc, model));
        update_expected_ss(var[i], &doc, ss);
        total += gsl_vector_get(likelihood, i);
	if ((verbose > 0) && ((i % (corpus->ndocs-1)) == 0) && (i>0)) 
	  Rprintf("lhood %5.5e   niter %5d\n", gsl_vector_get(likelihood, i), var[i]->niter);
        *avg_niter += var[i]->niter;
        *converged_pct += var[i]->converged;
        gsl_matrix_set_row(corpus_lambda, i, var[i]->lambda);
        gsl_matrix_set_row(corpus_nu, i, var[i]->nu);
        col_sum(var[i]->phi, phi_sum);
        gsl_matrix_set_row(corpus_phi_sum, i, phi_sum);
    }
    gsl_vector_free(phi_sum);
    *avg_niter = *avg_niter / corpus->ndocs;
    *converged_pct = *converged_pct / corpus->ndocs;
    *total_lhood = total;
}

/*
 * BG: copied from estimate.c
 * m step
 *
 */

void cov_shrinkage(gsl_matrix* mle, int n, gsl_matrix* result)
{
    int p = mle->size1, i;
    double temp = 0, alpha = 0, tau = 0, log_lambda_s = 0;
    gsl_vector
        *lambda_star = gsl_vector_calloc(p),
        t, u,
        *eigen_vals = gsl_vector_calloc(p),
        *s_eigen_vals = gsl_vector_calloc(p);
    gsl_matrix
        *d = gsl_matrix_calloc(p,p),
        *eigen_vects = gsl_matrix_calloc(p,p),
        *s_eigen_vects = gsl_matrix_calloc(p,p),
        *result1 = gsl_matrix_calloc(p,p);

    // get eigen decomposition

    sym_eigen(mle, eigen_vals, eigen_vects);
    for (i = 0; i < p; i++)
    {

        // compute shrunken eigenvalues

        temp = 0;
        alpha = 1.0/(n+p+1-2*i);
        vset(lambda_star, i, n * alpha * vget(eigen_vals, i));
    }

    // get diagonal mle and eigen decomposition

    t = gsl_matrix_diagonal(d).vector;
    u = gsl_matrix_diagonal(mle).vector;
    gsl_vector_memcpy(&t, &u);
    sym_eigen(d, s_eigen_vals, s_eigen_vects);

    // compute tau^2

    for (i = 0; i < p; i++)
        log_lambda_s += log(vget(s_eigen_vals, i));
    log_lambda_s = log_lambda_s/p;
    for (i = 0; i < p; i++)
        tau += pow(log(vget(lambda_star, i)) - log_lambda_s, 2)/(p + 4) - 2.0 / n;

    // shrink \lambda* towards the structured eigenvalues

    for (i = 0; i < p; i++)
        vset(lambda_star, i,
             exp((2.0/n)/((2.0/n) + tau) * log_lambda_s +
                 tau/((2.0/n) + tau) * log(vget(lambda_star, i))));

    // put the eigenvalues in a diagonal matrix

    t = gsl_matrix_diagonal(d).vector;
    gsl_vector_memcpy(&t, lambda_star);

    // reconstruct the covariance matrix

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, d, eigen_vects, 0, result1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, eigen_vects, result1, 0, result);

    // clean up

    gsl_vector_free(lambda_star);
    gsl_vector_free(eigen_vals);
    gsl_vector_free(s_eigen_vals);
    gsl_matrix_free(d);
    gsl_matrix_free(eigen_vects);
    gsl_matrix_free(s_eigen_vects);
    gsl_matrix_free(result1);
}



void maximization(llna_model* model, llna_ss* ss)
{
    int i, j;
    double sum;

    // mean maximization

    for (i = 0; i < model->k-1; i++)
        vset(model->mu, i, vget(ss->mu_ss, i) / ss->ndata);

    // covariance maximization

    for (i = 0; i < model->k-1; i++)
    {
        for (j = 0; j < model->k-1; j++)
        {
            mset(model->cov, i, j,
                 (1.0 / ss->ndata) *
                 (mget(ss->cov_ss, i, j) +
                  ss->ndata * vget(model->mu, i) * vget(model->mu, j) -
                  vget(ss->mu_ss, i) * vget(model->mu, j) -
                  vget(ss->mu_ss, j) * vget(model->mu, i)));
        }
    }
    if (PARAMS.cov_estimate == 1)
    {
        cov_shrinkage(model->cov, ss->ndata, model->cov);
    }
    matrix_inverse(model->cov, model->inv_cov);
    model->log_det_inv_cov = log_det(model->inv_cov);

    // topic maximization

    for (i = 0; i < model->k; i++)
    {
        sum = 0;
        for (j = 0; j < model->log_beta->size2; j++)
            sum += mget(ss->beta_ss, i, j);

        if (sum == 0) sum = safe_log(sum) * model->log_beta->size2;
        else sum = safe_log(sum);

        for (j = 0; j < model->log_beta->size2; j++)
            mset(model->log_beta, i, j, safe_log(mget(ss->beta_ss, i, j)) - sum);
    }
}

corpus* DocumentTermMatrix2Corpus(int *i, int *j, double *v, int nrow, int ncol, int length)
{
  corpus* c;
  int l, k;
  int* words;
  
  c = malloc(sizeof(corpus));
  c->nterms = ncol;
  c->ndocs = nrow;
  c->docs = (doc*) malloc(sizeof(doc)*(c->ndocs));
  words = malloc(sizeof(int)*nrow);  
  for (l = 0; l < nrow; l++) {
    c->docs[l].total = 0;
    c->docs[l].nterms = 0;
    words[l] = 0;
  }
  for (k = 0; k < length; k++) {
    c->docs[i[k]-1].nterms += 1;
  }
  for (l = 0; l < c->ndocs; l++) {
    c->docs[l].word = malloc(sizeof(int)*(c->docs[l].nterms));
    c->docs[l].count = malloc(sizeof(int)*(c->docs[l].nterms));
  }
  for (k = 0; k < length; k++) {
    c->docs[i[k]-1].word[words[i[k]-1]] = j[k] - 1;
    c->docs[i[k]-1].count[words[i[k]-1]] = v[k];
    c->docs[i[k]-1].total += v[k];
    words[i[k]-1] += 1;
  }
  return(c);
}

/*
 * BG: copied from estimate.c
 * run em
 *
 */

llna_model* em_initial_model(int k, corpus* corpus, const char* start, SEXP init_model)
{
    llna_model* model;
    if (PARAMS.verbose > 0) printf("starting from %s\n", start);
    if (strcmp(start, "rand")==0)
      model = random_init(k, corpus->nterms, PARAMS.verbose, PARAMS.seed);
    else if (strcmp(start, "seed")==0)
      model = corpus_init(k, corpus, PARAMS.verbose, PARAMS.seed);
    else
        model = R2llna_model(init_model);
    return(model);
}

SEXP rctm(SEXP i, SEXP j, SEXP v, SEXP nrow, SEXP ncol,
	  SEXP control, SEXP k, SEXP prefix, SEXP init_model) 
{
  corpus* corpus;
  llna_model *model;
  llna_ss* ss;
  FILE* lhood_fptr = NULL;
  char string[100];
  double convergence = 1, lhood = 0, lhood_old = 0;
  time_t t1,t2;
  double avg_niter, converged_pct, old_conv = 0;
  gsl_vector *likelihood;
  double *logLiks;
  gsl_matrix *corpus_lambda, *corpus_nu, *corpus_phi_sum;
  short reset_var = 1;
  const char *start, *dir;
  int NTOPICS, d, iteration, keep_iter, verbose;
  llna_var_param **var;
  SEXP ans;
     
  PARAMS.var_max_iter = *INTEGER(GET_SLOT(GET_SLOT(control, install("var")), install("iter.max")));
  PARAMS.var_convergence = *REAL(GET_SLOT(GET_SLOT(control, install("var")), install("tol")));
  PARAMS.em_max_iter = *INTEGER(GET_SLOT(GET_SLOT(control, install("em")), install("iter.max")));
  PARAMS.em_convergence = *REAL(GET_SLOT(GET_SLOT(control, install("em")), install("tol")));
  PARAMS.cg_max_iter = *INTEGER(GET_SLOT(GET_SLOT(control, install("cg")), install("iter.max")));
  PARAMS.cg_convergence = *REAL(GET_SLOT(GET_SLOT(control, install("cg")), install("tol")));
  PARAMS.verbose = *INTEGER(GET_SLOT(control, install("verbose")));
  PARAMS.save = *INTEGER(GET_SLOT(control, install("save")));
  PARAMS.seed = *INTEGER(GET_SLOT(control, install("seed")));
  PARAMS.keep = *INTEGER(GET_SLOT(control, install("keep")));
  /*   No shrinkage estimation is made because it is unclear what this option does. 
       The estimator based on a simple hierarchical model proposed in
       Daniels and Kass (2001) is used.  But sometimes rescaled eigenvalues
       are used where weights depend on the index but ordering of the
       eigenvalues is not ensured. */
  //  PARAMS.cov_estimate = *LOGICAL(GET_SLOT(control, install("shrinkage.covariance")));//
  PARAMS.cov_estimate = 0;

  NTOPICS = *INTEGER(AS_INTEGER(k));

  corpus = DocumentTermMatrix2Corpus(INTEGER(i),
				     INTEGER(j),
				     REAL(v),
				     *INTEGER(nrow),
				     *INTEGER(ncol),
				     LENGTH(v));
  start = CHAR(asChar(GET_SLOT(control, install("initialize"))));
  dir = CHAR(asChar(prefix));

  var = malloc(sizeof(llna_var_param*) * (corpus->ndocs));
  likelihood = gsl_vector_alloc(corpus->ndocs); 
  // START CODE FROM em (estimate.c)
  // set up the log likelihood log file
  
  if (PARAMS.save > 0) {
    sprintf(string, "%s/likelihood.dat", dir);
    lhood_fptr = fopen(string, "w");
  }
  if ((PARAMS.keep > 0) && (PARAMS.em_max_iter > 0)) {
    logLiks = malloc(sizeof(double*)*(ceil(PARAMS.em_max_iter/PARAMS.keep)));
  }

  // run em

  model = em_initial_model(NTOPICS, corpus, start, init_model);
  ss = new_llna_ss(model);
  corpus_lambda = gsl_matrix_alloc(corpus->ndocs, model->k);
  corpus_nu = gsl_matrix_alloc(corpus->ndocs, model->k);
  corpus_phi_sum = gsl_matrix_alloc(corpus->ndocs, model->k);
  (void) time(&t1);
  init_temp_vectors(model->k-1); // !!! hacky
  keep_iter = iteration = 0;
  sprintf(string, "%s/start", dir);
  if (PARAMS.save > 0) write_llna_model(model, string, PARAMS.verbose);
  while ((iteration < PARAMS.em_max_iter) &&
	 ((convergence > PARAMS.em_convergence) || (convergence < 0)))
    {
      verbose = PARAMS.verbose > 0 && (iteration % PARAMS.verbose) == 0;
      if (verbose) Rprintf("***** EM ITERATION %d *****\n", iteration+1);
      
      expectation(corpus, model, ss, &avg_niter, &lhood,
		  likelihood, 
		  corpus_lambda, corpus_nu, corpus_phi_sum,
		  reset_var, &converged_pct, var,
		  verbose);
      convergence = (lhood_old - lhood) / lhood_old;

      if (PARAMS.save > 0) {
	(void) time(&t2);
	fprintf(lhood_fptr, "%d %5.5e %5.5e %5.0f %5.5f %1.5f\n",
		iteration, lhood, convergence, difftime(t2, t1), avg_niter, converged_pct);
	if (((iteration % PARAMS.save)==0) || isnan(lhood))
	  {
	    sprintf(string, "%s/%03d", dir, iteration);
	    write_llna_model(model, string, PARAMS.verbose);
	    sprintf(string, "%s/%03d-lambda.dat", dir, iteration);
	    printf_matrix(string, corpus_lambda);
	    sprintf(string, "%s/%03d-nu.dat", dir, iteration);
	    printf_matrix(string, corpus_nu);
	  }
	(void) time(&t1);
      }
      if (convergence < 0)
        {
	  reset_var = 0;
	  if (PARAMS.var_max_iter > 0)
	    PARAMS.var_max_iter += 10;
	  else PARAMS.var_convergence /= 10;
        }
      else
        {
	  if ((PARAMS.keep > 0) && (PARAMS.em_max_iter > 0)) {
	    logLiks[keep_iter] = lhood;
	    keep_iter++;
	  }
	  maximization(model, ss);
	  lhood_old = lhood;
	  reset_var = 1;
	  iteration++;
        }
      
      if (PARAMS.save > 0) fflush(lhood_fptr);
      reset_llna_ss(ss);
      old_conv = convergence;
    }

  if (PARAMS.em_max_iter < 0) {
    expectation(corpus, model, ss, &avg_niter, &lhood,
		likelihood, 
		corpus_lambda, corpus_nu, corpus_phi_sum,
		reset_var, &converged_pct, var,
		verbose);
  }


  if (PARAMS.save > 0) {
    sprintf(string, "%s/final", dir);
    write_llna_model(model, string, PARAMS.verbose);
    sprintf(string, "%s/final-lambda.dat", dir);
    printf_matrix(string, corpus_lambda);
    sprintf(string, "%s/final-nu.dat", dir);
    printf_matrix(string, corpus_nu);
    fclose(lhood_fptr);
  }
  // END CODE FROM em (estimate.c)
  if ((PARAMS.keep > 0) && (PARAMS.em_max_iter > 0)) {
    logLiks = realloc(logLiks, sizeof(double*)*(keep_iter));
  }

  // construct return object
  PROTECT(ans = NEW_OBJECT(MAKE_CLASS("CTM_VEM")));
  ans = returnObjectCTM(ans, model, corpus, corpus_lambda, corpus_nu, corpus_phi_sum, var, likelihood, iteration, logLiks, keep_iter);

  for (d = 0; d < corpus->ndocs; d++) free_llna_var_param(var[d]);
  free(var); free(corpus); free(ss); free(model);
  gsl_matrix_free(corpus_lambda); gsl_matrix_free(corpus_nu); gsl_matrix_free(corpus_phi_sum);
  gsl_vector_free(likelihood);
  UNPROTECT(1);
  return(ans);
}
