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

#include <stdio.h>
#include <stdlib.h>
#include "dataset.h"

using namespace std;

int dataset::readDocumentTermMatrix(int *i, int *j, double *v, int length) {

  int k, l, m;
  int *words, *lengths;
  document **documents;

  words = new int[M];
  lengths = new int[M];
  documents = new document*[M];
  // initialize lengths
  for (l = 0; l < M; l++) {
    lengths[l] = 0;
    words[l] = 0;
  }
  // determine total number of words in each document
  for (k = 0; k < length; k++) {
    lengths[i[k]-1] += v[k];
  } 
  // allocate words
  for (l = 0; l < M; l++) {
    documents[l] = new document(lengths[l]);
  } 
  // assign words in documents
  for (k = 0; k < length; k++) {
    for (m = 0; m < v[k]; m++) {
      documents[i[k]-1]->words[words[i[k]-1]] = j[k] - 1;
      words[i[k]-1] += 1;
    }
  }
  for (l = 0; l < M; l++) {
    add_doc(documents[l], l);
  } 
  return 0;
}
