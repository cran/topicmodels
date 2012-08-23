// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)
// (C) Copyright 2009, Bettina Gruen (Bettina [dot] Gruen [at] wu-wien [dot] ac [dot] at)

// This file re-uses parts of LDA-C, especially of the function run_em.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "rlda.h"

/*
 * ArgMax
 *
 */

int ArgMax(double* x, int n)
{
    int i;
    double max = x[0];
    int ArgMax = 0;
    for (i = 1; i < n; i++)
    {
        if (x[i] > max)
        {
            max = x[i];
            ArgMax = i;
        }
    }
    return(ArgMax);
}

/*
 * BG: max_corpus_length copied from lda-data.c
 *
 */

int max_corpus_length(corpus* c)
{
    int n, max = 0;
    for (n = 0; n < c->num_docs; n++)
	if (c->docs[n].length > max) max = c->docs[n].length;
    return(max);
}

/*
 * BG: save_gamma copied from lda-estimate.c
 * saves the gamma parameters of the current dataset
 *
 */

void save_gamma(char* filename, double** gamma, int num_docs, int num_topics)
{
    FILE* fileptr;
    int d, k;
    fileptr = fopen(filename, "w");

    for (d = 0; d < num_docs; d++)
    {
	fprintf(fileptr, "%5.10f", gamma[d][0]);
	for (k = 1; k < num_topics; k++)
	{
	    fprintf(fileptr, " %5.10f", gamma[d][k]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}


/*
 * BG: doc_e_step copied from lda-estimate.c
 * perform inference on a document and update sufficient statistics
 *
 */

double doc_e_step(document* doc, double* gamma, double** phi,
                  lda_model* model, lda_suffstats* ss)
{
    double likelihood;
    int n, k;

    // posterior inference

    likelihood = lda_inference(doc, model, gamma, phi);

    // update sufficient statistics

    double gamma_sum = 0;
    for (k = 0; k < model->num_topics; k++)
    {
        gamma_sum += gamma[k];
        ss->alpha_suffstats += digamma(gamma[k]);
    }
    ss->alpha_suffstats -= model->num_topics * digamma(gamma_sum);

    for (n = 0; n < doc->length; n++)
    {
        for (k = 0; k < model->num_topics; k++)
        {
            ss->class_word[k][doc->words[n]] += doc->counts[n]*phi[n][k];
            ss->class_total[k] += doc->counts[n]*phi[n][k];
        }
    }

    ss->num_docs = ss->num_docs + 1;

    return(likelihood);
}

/*
 * BG: construct return object
 *
 */

SEXP returnObjectLDA(SEXP ans, lda_model* model, corpus* corpus, double ***phi, 
		     double **var_gamma, double *likelihood, int iter, double *logLiks, int keep_iter) {
  SEXP tp, I, J, V, wordassign, nms, dn;
  SEXP ROWNAMES, COLNAMES;
  document* doc;
  int i, j, d, total, n;
  double *m;

  tp = PROTECT(allocVector(INTSXP, 1));
  *INTEGER(tp) = iter;
  SET_SLOT(ans, install("iter"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocVector(INTSXP, 1));
  *INTEGER(tp) = model->num_topics;
  SET_SLOT(ans, install("k"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocVector(REALSXP, 1));
  *REAL(tp) = model->alpha;
  SET_SLOT(ans, install("alpha"), tp);
  UNPROTECT(1);
  
  tp = PROTECT(allocVector(REALSXP, corpus->num_docs));
  for (i = 0; i < corpus->num_docs; i++)
    REAL(tp)[i] = likelihood[i];
  SET_SLOT(ans, install("loglikelihood"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocVector(INTSXP, 2));
  INTEGER(tp)[0] = corpus->num_docs;
  INTEGER(tp)[1] = corpus->num_terms;
  SET_SLOT(ans, install("Dim"), tp);
  UNPROTECT(1);
  
  tp = PROTECT(allocMatrix(REALSXP, model->num_topics, corpus->num_terms));
  for (i = 0; i < model->num_topics; i++)
    for (j = 0; j < corpus->num_terms; j++)
      REAL(tp)[i + model->num_topics * j] = model->log_prob_w[i][j];
  SET_SLOT(ans, install("beta"), tp);
  UNPROTECT(1);

  tp = PROTECT(allocMatrix(REALSXP, corpus->num_docs, model->num_topics));
  m = REAL(tp);
  for (i = 0; i < corpus->num_docs; i++)
    for (j = 0; j < model->num_topics; j++)
      m[i + corpus->num_docs * j] = var_gamma[i][j];
  SET_SLOT(ans, install("gamma"), tp);
  UNPROTECT(1);
   
  if ((KEEP > 0) && (EM_MAX_ITER > 0)) {
    tp = PROTECT(allocVector(REALSXP, keep_iter));
    for (i = 0; i < keep_iter; i++) 
      REAL(tp)[i] = logLiks[i];
    SET_SLOT(ans, install("logLiks"), tp);
    UNPROTECT(1);
  }
  
  wordassign = PROTECT(allocVector(VECSXP, 6));
  total = 0;
  for (d = 0; d < corpus->num_docs; d++) {
    doc = &(corpus->docs[d]);
    total += doc->length;
  }

  I = PROTECT(allocVector(INTSXP, total));
  J = PROTECT(allocVector(INTSXP, total));
  V = PROTECT(allocVector(REALSXP, total));

  i = 0;
  for (d = 0; d < corpus->num_docs; d++) {
    doc = &(corpus->docs[d]);
    for (n = 0; n < doc->length; n++) {
      INTEGER(I)[i] = d+1;
      INTEGER(J)[i] = doc->words[n] + 1;
      REAL(V)[i] = ArgMax(phi[d][n], model->num_topics) + 1;
      i++;
    }
  }
  SET_VECTOR_ELT(wordassign, 0, I); 
  SET_VECTOR_ELT(wordassign, 1, J);
  SET_VECTOR_ELT(wordassign, 2, V);
  UNPROTECT(3);

  tp = PROTECT(allocVector(INTSXP, 1));
  INTEGER(tp)[0] = corpus->num_docs;
  SET_VECTOR_ELT(wordassign, 3, tp);
  UNPROTECT(1);
  tp = PROTECT(allocVector(INTSXP, 1));
  INTEGER(tp)[0] = corpus->num_terms;
  SET_VECTOR_ELT(wordassign, 4, tp);
  UNPROTECT(1);

  dn = PROTECT(allocVector(VECSXP, 2));
  ROWNAMES = PROTECT(allocVector(INTSXP, corpus->num_docs));
  COLNAMES = PROTECT(allocVector(INTSXP, corpus->num_terms));
  for (d = 0; d < corpus->num_docs; d++) INTEGER(ROWNAMES)[d] = d+1;
  for (d = 0; d < corpus->num_terms; d++) INTEGER(COLNAMES)[d] = d+1;
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

lda_model* R2lda_model(SEXP init_model)
{
    int i, j, num_terms, num_topics;
    lda_model* model;

    num_terms = INTEGER(GET_SLOT(init_model, install("Dim")))[1];
    num_topics = *INTEGER(GET_SLOT(init_model, install("k")));
    model = new_lda_model(num_terms, num_topics);
    model->alpha = *REAL(GET_SLOT(init_model, install("alpha")));
     
    for (i = 0; i < num_topics; i++)
      for (j = 0; j < num_terms; j++) 
	model->log_prob_w[i][j] = REAL(GET_SLOT(init_model, install("beta")))[i + num_topics * j];
    return(model);
}

corpus* DocumentTermMatrix2corpus(int *i, int *j, double *v, int nrow, int ncol, int length)
{
  corpus* c;
  int l, k;
  int* words;
  
  c = malloc(sizeof(corpus));
  c->num_terms = ncol;
  c->num_docs = nrow;
  words = malloc(sizeof(int)*nrow);
  c->docs = (document*) malloc(sizeof(document)*nrow);
  for (l = 0; l < nrow; l++) {
    c->docs[l].total = 0;
    c->docs[l].length = 0;
    words[l] = 0;
  }
  for (k = 0; k < length; k++) {
    c->docs[i[k]-1].length += 1;
  }
  for (l = 0; l < c->num_docs; l++) {
    c->docs[l].words = malloc(sizeof(int)*(c->docs[l].length));
    c->docs[l].counts = malloc(sizeof(int)*(c->docs[l].length));
  }
  for (k = 0; k < length; k++) {
    c->docs[i[k]-1].words[words[i[k]-1]] = j[k] - 1;
    c->docs[i[k]-1].counts[words[i[k]-1]] = v[k];
    c->docs[i[k]-1].total += v[k];
    words[i[k]-1] += 1;
  }
  free(words);
  return(c);
}

SEXP rlda(SEXP i, SEXP j, SEXP v, SEXP nrow, SEXP ncol,
          SEXP control, SEXP k, SEXP prefix, SEXP init_model) 
{

  corpus* corpus;
  const char *start, *directory;
  SEXP ans; 
  int iter, keep_iter, d, n, max_length, verbose = 0;
  lda_model *model = NULL;
  double **var_gamma, ***phi, *llh, *logLiks = NULL;
  FILE* likelihood_file = NULL;
  char filename[100];
  lda_suffstats* ss = NULL;
  double likelihood, likelihood_old = 0, converged = 1;

  VAR_MAX_ITER = *INTEGER(GET_SLOT(GET_SLOT(control, install("var")), install("iter.max")));
  VAR_CONVERGED = *REAL(GET_SLOT(GET_SLOT(control, install("var")), install("tol")));
  EM_MAX_ITER = *INTEGER(GET_SLOT(GET_SLOT(control, install("em")), install("iter.max")));
  EM_CONVERGED = *REAL(GET_SLOT(GET_SLOT(control, install("em")), install("tol")));
  LAG = *INTEGER(GET_SLOT(control, install("verbose")));
  SAVE = *INTEGER(GET_SLOT(control, install("save")));
  KEEP = *INTEGER(GET_SLOT(control, install("keep")));
  ESTIMATE_ALPHA = *LOGICAL(GET_SLOT(control, install("estimate.alpha")));
  SEED = *INTEGER(GET_SLOT(control, install("seed")));
  seedMT(SEED);

  NTOPICS = *INTEGER(AS_INTEGER(k));
  INITIAL_ALPHA = *REAL(GET_SLOT(control, install("alpha")));
  corpus = DocumentTermMatrix2corpus(INTEGER(i),
				     INTEGER(j),
				     REAL(v),
				     *INTEGER(nrow),
				     *INTEGER(ncol),
				     LENGTH(v));
  start = CHAR(asChar(GET_SLOT(control, install("initialize"))));
  directory = CHAR(asChar(prefix));
  
  llh = malloc(sizeof(double)*(corpus->num_docs));
  // START CODE FROM run_em (lda-estimate.c)
  // allocate variational parameters
  // BG: allocate for phi[d][n][i] 

  max_length = max_corpus_length(corpus);
  var_gamma = malloc(sizeof(double*)*(corpus->num_docs));
  phi = malloc(sizeof(double**)*(corpus->num_docs));
  if ((KEEP > 0) && (EM_MAX_ITER > 0)) {
    logLiks = malloc(sizeof(double*)*(ceil((double)EM_MAX_ITER/KEEP)));
  }
  for (d = 0; d < corpus->num_docs; d++) {
    var_gamma[d] = malloc(sizeof(double) * NTOPICS);
    phi[d] = malloc(sizeof(double*)*max_length);
    for (n = 0; n < max_length; n++)
      phi[d][n] = malloc(sizeof(double) * NTOPICS);
  }

  // initialize model
  
  if (strcmp(start, "seeded")==0)
    {
      model = new_lda_model(corpus->num_terms, NTOPICS);
      ss = new_lda_suffstats(model);
      corpus_initialize_ss(ss, model, corpus);
      lda_mle(model, ss, 0, LAG > 0);
      model->alpha = INITIAL_ALPHA;
    }
  else if (strcmp(start, "random")==0)
    {
      model = new_lda_model(corpus->num_terms, NTOPICS);
      ss = new_lda_suffstats(model);
      random_initialize_ss(ss, model);
      lda_mle(model, ss, 0, LAG > 0);
      model->alpha = INITIAL_ALPHA;
    }
  else 
    {
      model = R2lda_model(init_model);
      ss = new_lda_suffstats(model);
    }
  
  sprintf(filename,"%s/000",directory);
  if (SAVE > 0) save_lda_model(model, filename);
  
  // run expectation maximization
  
  iter = 0;
  keep_iter = 0;

  if (SAVE > 0) {
    sprintf(filename, "%s/likelihood.dat", directory);
    likelihood_file = fopen(filename, "w");
  }
  
  while (((converged < 0) || (converged > EM_CONVERGED) || (iter <= 2)) && (iter < EM_MAX_ITER))
    {
      iter++; 
      verbose = (LAG > 0) && ((iter % LAG) == 0);
      if (verbose) Rprintf("**** em iteration %d ****\n", iter);
      likelihood = 0;
      zero_initialize_ss(ss, model);
      
      // e-step
      
      for (d = 0; d < corpus->num_docs; d++)
        {
	  if (verbose && ((d % 1000) == 0) && (d>0)) Rprintf("document %d\n",d);
	  likelihood += doc_e_step(&(corpus->docs[d]),
				   var_gamma[d],
				   phi[d],
				   model,
				   ss);
        }
      
      // m-step
      
      lda_mle(model, ss, ESTIMATE_ALPHA, verbose);
      
      // check for convergence
      
      converged = (likelihood_old - likelihood) / (likelihood_old);
      if (converged < 0) VAR_MAX_ITER = VAR_MAX_ITER * 2;
      likelihood_old = likelihood;
      
      if ((KEEP > 0) && (EM_MAX_ITER > 0)) {
	logLiks[keep_iter] = likelihood;
	keep_iter++;
      }
		
      // output model and likelihood
      
      if (SAVE > 0) {
	fprintf(likelihood_file, "%10.10f\t%5.5e\n", likelihood, converged);
	fflush(likelihood_file);
	if ((iter % SAVE) == 0)
	  {
	    sprintf(filename,"%s/%03d",directory, iter);
	    save_lda_model(model, filename);
	    sprintf(filename,"%s/%03d.gamma",directory, iter);
	    save_gamma(filename, var_gamma, corpus->num_docs, model->num_topics);
	  }
      }
    }

  // output the word assignments (for visualization)
  
  for (d = 0; d < corpus->num_docs; d++) {
      if (verbose && ((d % 100) == 0) && (d>0)) Rprintf("final e step document %d\n",d);
      llh[d] = lda_inference(&(corpus->docs[d]), model, var_gamma[d], phi[d]);
      likelihood += llh[d];
    }
  // END CODE FROM run_em (lda-estimate.c)
  // BG: word assignments are not saved but returned in the R object
  // output the final model
  
  if (SAVE > 0) {
    fclose(likelihood_file);
    sprintf(filename,"%s/final",directory);
    save_lda_model(model, filename);
    sprintf(filename,"%s/final.gamma",directory);
    save_gamma(filename, var_gamma, corpus->num_docs, model->num_topics);
  }
  
  if ((KEEP > 0) && (EM_MAX_ITER > 0)) {
    logLiks = realloc(logLiks, sizeof(double*)*(keep_iter));
  }

  // construct return object
  PROTECT(ans = NEW_OBJECT(MAKE_CLASS("LDA_VEM")));
  ans = returnObjectLDA(ans, model, corpus, phi, var_gamma, llh, iter, logLiks, keep_iter);

  for (d = 0; d < corpus->num_docs; d++) {
    free(var_gamma[d]);
    for (n = 0; n < max_length; n++) {
      free(phi[d][n]);
    }
    free(phi[d]);
  }
  free(phi); free(var_gamma); free(llh); free(logLiks);
  free(corpus); free_lda_suffstats(ss, model->num_topics, model->num_terms); free_lda_model(model); 
  UNPROTECT(1);
  return(ans);
}
