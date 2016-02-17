/*
 * Original code by Xuan-Hieu Phan modified by Bettina Gruen
 *
 * Copyright (C) 2007 by
 * 
 * 	Xuan-Hieu Phan
 *	hieuxuan@ecei.tohoku.ac.jp or pxhieu@gmail.com
 * 	Graduate School of Information Sciences
 * 	Tohoku University
 *
 * Copyright (C) 2009 by
 * 
 * 	Bettina Gruen
 *
 * GibbsLDA++ is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * GibbsLDA++ is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/* 
 * References:
 * + The Java code of Gregor Heinrich (gregor@arbylon.net)
 *   http://www.arbylon.net/projects/LdaGibbsSampler.java
 * + "Parameter estimation for text analysis" by Gregor Heinrich
 *   http://www.arbylon.net/publications/text-est.pdf
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constants.h"
#include "utilities.h"
#include "dataset.h"
#include "model.h"

using namespace std;

model::~model() {
  int m, k, w;

    if (logLiks) {
      delete[] logLiks;
    }

    if (p) {
	delete[] p;
    }

    if (ptrndata) {
        delete ptrndata;
    }
    
    if (z) {
	for (m = 0; m < M; m++) {
	    if (z[m]) {
		delete[] z[m];
	    }
	}
	delete[] z;
    }
    if (wordassign) {
	for (m = 0; m < M; m++) {
	    if (wordassign[m]) {
		delete[] wordassign[m];
	    }
	}
	delete[] wordassign;
    }
    if (beta) {
      for (w = 0; w < V; w++) {
	if (beta[w]) {
	  delete[] beta[w];
	}
      }
      delete[] beta;
    }

    if (nw) {
      for (w = 0; w < V; w++) {
	if (nw[w]) {
	  delete[] nw[w];
	}
      }
      delete[] nw;
    }

    if (nd) {
	for (m = 0; m < M; m++) {
	    if (nd[m]) {
		delete[] nd[m];
	    }
	}
	delete[] nd;
    } 
    
    if (Vbeta) {
      delete[] Vbeta;
    }
    
    if (nwsum) {
	delete[] nwsum;
    }   
    if (ndsum) {
	delete[] ndsum;
    }
    if (theta) {
	for (m = 0; m < M; m++) {
	    if (theta[m]) {
		delete[] theta[m];
	    }
	}
	delete[] theta;
    }
    
    if (phi) {
	for (k = 0; k < K; k++) {
	    if (phi[k]) {
		delete[] phi[k];
	    }
	}
    }
    delete[] phi;
}

void model::set_default_values() {
    tassign_suffix = ".tassign";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    others_suffix = ".others";
    
    dir = "./";
    model_name = "model-final";    
    
    ptrndata = NULL;
    
    M = 0;
    V = 0;
    K = 100;
    alpha = 50.0 / K;
    beta1 = 0.1;
    niters = 2000;
    liter = 0;
    verbose = 200;    
    save = 0;
    keep = 0;
    loglikelihood = 0;
    estimate_phi = 1;

    beta = NULL;
    logLiks = NULL;
    p = NULL;
    z = NULL;
    wordassign = NULL;
    nw = NULL;
    nd = NULL;
    nwsum = NULL;
    ndsum = NULL;
    Vbeta = NULL;
    theta = NULL;
    phi = NULL;
}

int model::save_model(string new_model_name) {
  if (save_model_tassign(dir + "/" + new_model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_model_others(dir + "/" + new_model_name + others_suffix)) {
	return 1;
    }
    
    if (save_model_theta(dir + "/" + new_model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_model_phi(dir + "/" + new_model_name + phi_suffix)) {
	return 1;
    }
    
    return 0;
}

int model::save_model_tassign(string filename) {
    int i, j;
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        Rprintf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    // write docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {    
	for (j = 0; j < ptrndata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", ptrndata->docs[i]->words[j], z[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);

    return 0;
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        Rprintf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < M; i++) {
	for (int j = 0; j < K; j++) {
	    fprintf(fout, "%f ", theta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_phi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	Rprintf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < V; j++) {
	    fprintf(fout, "%f ", phi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	Rprintf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::init(int *i, int *j, int *v, int total, double *delta, int *Z, double *Phi, int init) {
    int m, n, w, k;
    
    if (verbose > 0) Rprintf("K = %d; V = %d; M = %d\n", K, V, M);

    p = new double[K];
    if (keep > 0) {
      int keep_iter = ceil((double)niters/keep);
      logLiks = new double[keep_iter];
    }

    // + read training data
    ptrndata = new dataset(M, V);
    ptrndata->readDocumentTermMatrix(i, j, v, total);
		
    // + allocate memory and assign values for variables
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, verbose: from command line or default values

    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }

    if (seeded == 1) {
      beta = new double*[V];
      for (w = 0; w < V; w++) {
        beta[w] = new double[K];
        for (k = 0; k < K; k++) {
	  beta[w][k] = delta[k + K * w];
        }
      }
      
      Vbeta = new double[K];
      for (int k = 0; k < K; k++) {
	Vbeta[k] = 0.0;
	for (w = 0; w < V; w++) {
	  Vbeta[k] += beta[w][k];
	}
      }
    } else {
      beta1 = delta[0]; 
    }

    z = new int*[M];
    wordassign = new int*[M];
    k = 0;
    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;
	z[m] = new int[N];
	wordassign[m] = new int[N];
	
        // initialize for z
        for (n = 0; n < N; n++) {
	  int topic = get_z(m, n, Phi, Z, k, init);
	  k++;
    	    z[m][n] = topic;
    	    // number of instances of word i assigned to topic j
    	    nw[ptrndata->docs[m]->words[n]][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }

    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    
    if (estimate_phi == 0) {
      for (int k = 0; k < K; k++) {
	for (int w = 0; w < V; w++) {
	  phi[k][w] = exp(Phi[k + K * w]);
	}
      }
    }
    return 0;
}

int model::get_z(int m, int n, double *Phi, int *Z, int index, int init) 
{
  int topic = 0;
  if (init == 0) {
    topic = (int)(K * unif_rand());
  }
  if (init == 1) {
    int w = ptrndata->docs[m]->words[n];
    for (int k = 0; k < K; k++) {
      p[k] = exp(Phi[k + K * w]);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
      p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = (unif_rand()) * p[K - 1];
    
    for (int k = 0; k < K; k++) {
      topic = k;
      if (p[k] > u) {
	break;
      }
    }
  }
  if (init == 2) {
    topic = Z[index] - 1;
  }
  return topic;
}

void model::estimate() {
  if (verbose > 0) Rprintf("Sampling %d iterations!\n", niters);
 
    int keep_iter = 0;
    int last_iter = liter;

    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
	// for all z_i
	for (int m = 0; m < M; m++) {
	    for (int n = 0; n < ptrndata->docs[m]->length; n++) {
		// (z_i = z[m][n])
		// sample from p(z_i|z_-i, w)
		z[m][n] = sampling(m, n);
	    }
	}
	
	if ((save > 0) && (liter % save == 0)) {
	  if (verbose > 0) Rprintf("Saving the model at iteration %d ...\n", liter);
	  compute_theta();
	  if (estimate_phi == 1) compute_phi();
	  save_model(utilities::generate_model_name(liter));
	  
	} else if ((verbose > 0) && (liter % verbose == 0)) Rprintf("Iteration %d ...\n", liter);
	if ((keep > 0) && (liter % keep == 0)) {
	  inference();
	  logLiks[keep_iter] = loglikelihood;
	  keep_iter++;
	}
    }
    
    if (verbose > 0) Rprintf("Gibbs sampling completed!\n");
    compute_theta();
    if (estimate_phi == 1) compute_phi();
    // compute wordassign;
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < ptrndata->docs[m]->length; n++) {
	wordassign[m][n] = get_wordassign(m, n);
      }
    }
    liter--;    
    if (save > 0) save_model(utilities::generate_model_name(-1));
}

int model::sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    nw[w][topic] -= 1;
    nd[m][topic] -= 1;
    nwsum[topic] -= 1;
    ndsum[m] -= 1;

    double Vbeta1 = V * beta1;
    double Kalpha = K * alpha;    

    if (estimate_phi == 1) {
      // do multinomial sampling via cumulative method
      if (seeded == 1) {
	for (int k = 0; k < K; k++) {
	  p[k] = (nw[w][k] + beta[w][k]) / (nwsum[k] + Vbeta[k]) *
	    (nd[m][k] + alpha) / (ndsum[m] + Kalpha);
	}
      } else {
	for (int k = 0; k < K; k++) {
	  p[k] = (nw[w][k] + beta1) / (nwsum[k] + Vbeta1) *
	    (nd[m][k] + alpha) / (ndsum[m] + Kalpha);
	}
      }
    } else {
      for (int k = 0; k < K; k++) {
	p[k] = phi[k][w] *
	  (nd[m][k] + alpha) / (ndsum[m] + Kalpha);
      }
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = unif_rand() * p[K - 1];
    
    for (int k = 0; k < K; k++) {
        topic = k;
	if (p[topic] > u) {
	    break;
	}
    }
    
    // add newly estimated z_i to count variables
    nw[w][topic] += 1;
    nd[m][topic] += 1;
    nwsum[topic] += 1;
    ndsum[m] += 1;    
    
    return topic;
}

int model::get_wordassign(int m, int n) {
    int topic = 0;
    int k, w;

    w = ptrndata->docs[m]->words[n];
    double Kalpha = K * alpha;    

    for (k = 0; k < K; k++) {
      p[k] = phi[k][w] *
	(nd[m][k] + alpha) / (ndsum[m] + Kalpha);
    }
    double pmax = 0.0;
    for (k = 0; k < K; k++) {
        if (p[k] > pmax) {
	  pmax = p[k];
	  topic = k;
	}
    }
    return topic;
}

void model::compute_theta() {
    for (int m = 0; m < M; m++) {
	for (int k = 0; k < K; k++) {
	    theta[m][k] = (nd[m][k] + alpha) / (ndsum[m] + K * alpha);
	}
    }
}

void model::compute_phi() {
  if (seeded == 1) {
    for (int k = 0; k < K; k++) {
      for (int w = 0; w < V; w++) {
	phi[k][w] = (nw[w][k] + beta[w][k]) / (nwsum[k] + Vbeta[k]);
      }
    }
  } else {
    double Vbeta1 = V * beta1;
    for (int k = 0; k < K; k++) {
      for (int w = 0; w < V; w++) {
	phi[k][w] = (nw[w][k] + beta1) / (nwsum[k] + Vbeta1);
      }
    }
  }
}

void model::inference() {
  if (seeded == 1) {
    loglikelihood = 0;
    for (int k = 0; k < K; k++) {
      loglikelihood += lgamma(Vbeta[k]);
      for (int w = 0; w < V; w++) {
	loglikelihood -= lgamma(beta[w][k]);
	loglikelihood += lgamma(nw[w][k] + beta[w][k]);
      }
      loglikelihood -= lgamma(nwsum[k] + Vbeta[k]);
    }
  } else {
    double Vbeta1 = V * beta1;
    loglikelihood = K * (lgamma(Vbeta1) - V * lgamma(beta1));
    for (int k = 0; k < K; k++) {
      for (int w = 0; w < V; w++) {
	loglikelihood += lgamma(nw[w][k] + beta1);
      }
      loglikelihood -= lgamma(nwsum[k] + Vbeta1);
    }
  }
}
