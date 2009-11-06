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
#include <string>
#include "utilities.h"
#include "model.h"

using namespace std;

string utilities::generate_model_name(int iter) {
    string model_name = "model-";

    char buff[BUFF_SIZE_SHORT];
    
    if (0 <= iter && iter < 10) {
	sprintf(buff, "0000%d", iter);
    } else if (10 <= iter && iter < 100) {
	sprintf(buff, "000%d", iter);
    } else if (100 <= iter && iter < 1000) {
	sprintf(buff, "00%d", iter);
    } else if (1000 <= iter && iter < 10000) {
	sprintf(buff, "0%d", iter);
    } else {
	sprintf(buff, "%d", iter);
    }
    
    if (iter >= 0) {
	model_name += buff;
    } else {
	model_name += "final";
    }
    
    return model_name;
}

