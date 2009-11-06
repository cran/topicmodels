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

#define random rand
#define srandom srand 

using namespace std;

model::~model() {
  int m, k, w;
    if (p) {
	delete[] p;
    }

//     if (ptrndata) {
//         delete ptrndata;
//     }
    
    if (z) {
	for (m = 0; m < M; m++) {
	    if (z[m]) {
		delete[] z[m];
	    }
	}
    }
    if (wordassign) {
	for (m = 0; m < M; m++) {
	    if (wordassign[m]) {
		delete[] wordassign[m];
	    }
	}
    }

    if (nd) {
	for (m = 0; m < M; m++) {
	    if (nd[m]) {
		delete[] nd[m];
	    }
	}
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
    }
    
    if (phi) {
	for (k = 0; k < K; k++) {
	    if (phi[k]) {
		delete[] phi[k];
	    }
	}
    }

    if (nw) {
      for (w = 0; w < V; w++) {
	if (nw[w]) {
	  delete[] nw[w];
	}
      }
    }
}

void model::set_default_values() {
    tassign_suffix = ".tassign";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    others_suffix = ".others";
    
    dir = "./";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    
    ptrndata = NULL;
    
    M = 0;
    V = 0;
    K = 100;
    alpha = 50.0 / K;
    beta = 0.1;
    niters = 2000;
    liter = 0;
    savestep = 200;    
    
    p = NULL;
    z = NULL;
    wordassign = NULL;
    nw = NULL;
    nd = NULL;
    nwsum = NULL;
    ndsum = NULL;
    theta = NULL;
    phi = NULL;
}

int model::save_model(string model_name) {
    if (save_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_model_theta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_model_phi(dir + model_name + phi_suffix)) {
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
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::init(int *i, int *j, double *v, int length) {
    int m, n, w, k;
    
    Rprintf("K= %d; V = %d; M = %d\n", K, V, M);
    p = new double[K];

    // + read training data
    ptrndata = new dataset(M, V);
    ptrndata->readDocumentTermMatrix(i, j, v, length);
		
    // + allocate memory and assign values for variables
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

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
    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    wordassign = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;
	z[m] = new int[N];
	wordassign[m] = new int[N];
	
        // initialize for z
        for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
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
    
    return 0;
}

int model::initc(int *i, int *j, double *v, int length, double *Phi) {
    int m, n, w, k;
    
    Rprintf("K= %d; V = %d; M = %d\n", K, V, M);
    p = new double[K];

    // + read training data
    ptrndata = new dataset(M, V);
    ptrndata->readDocumentTermMatrix(i, j, v, length);
		
    // + allocate memory and assign values for variables
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

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
    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    wordassign = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;
	z[m] = new int[N];
	wordassign[m] = new int[N];
	
        // initialize for z
        for (n = 0; n < N; n++) {
	    int topic = get_z(m, n, Phi);
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
    
    return 0;
}

int model::get_z(int m, int n, double *Phi) 
{
  int topic;
  int w = ptrndata->docs[m]->words[n];
  for (int k = 0; k < K; k++) {
    p[k] = Phi[k * (w  + 1)];
  }
  // cumulate multinomial parameters
  for (int k = 1; k < K; k++) {
    p[k] += p[k - 1];
  }
  // scaled sample because of unnormalized p[]
  double u = ((double)random() / RAND_MAX) * p[K - 1];
  
  for (int k = 0; k < K; k++) {
    topic = k;
    if (p[k] > u) {
      break;
    }
  }
  return topic;
}

void model::estimate() {
    Rprintf("Sampling %d iterations!\n", niters);

    int last_iter = liter;
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
        // Rprintf("Iteration %d ...\n", liter);
	
	// for all z_i
	for (int m = 0; m < M; m++) {
	    for (int n = 0; n < ptrndata->docs[m]->length; n++) {
		// (z_i = z[m][n])
		// sample from p(z_i|z_-i, w)
		z[m][n] = sampling(m, n);
	    }
	}
	
	if (savestep > 0) {
	    if (liter % savestep == 0) {
		// saving the model
		Rprintf("Saving the model at iteration %d ...\n", liter);
		compute_theta();
		compute_phi();
		save_model(utilities::generate_model_name(liter));
	    }
	}
    }
    
    Rprintf("Gibbs sampling completed!\n");
    compute_theta();
    compute_phi();
    // compute wordassign;
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < ptrndata->docs[m]->length; n++) {
	wordassign[m][n] = get_wordassign(m, n);
      }
    }
    liter--;    
    save_model(utilities::generate_model_name(-1));
}

int model::sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    nw[w][topic] -= 1;
    nd[m][topic] -= 1;
    nwsum[topic] -= 1;
    ndsum[m] -= 1;

    double Vbeta = V * beta;
    double Kalpha = K * alpha;    
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
	p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
		    (nd[m][k] + alpha) / (ndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
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
    int topic;
    int w;

    w = ptrndata->docs[m]->words[n];
    double Vbeta = V * beta;
    double Kalpha = K * alpha;    

    for (int k = 0; k < K; k++) {
	p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
		    (nd[m][k] + alpha) / (ndsum[m] + Kalpha);
    }
    double pmax = 0.0;
    for (int k = 0; k < K; k++) {
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
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < V; w++) {
	    phi[k][w] = (nw[w][k] + beta) / (nwsum[k] + V * beta);
	}
    }
}
