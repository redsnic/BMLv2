/*
 * SearchDAGs.h
 *
 *  Created on: 24 ago 2017
 *      Author: rossi
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>

#ifndef SEARCHDAGS_H_

#ifndef DS
#include "DataStructures.h"
#define DS
#endif

#include "Utils.h"

void dfs(vector<bool>& conf, vector<bool>& seen, int& label, vector< vector< double> >& cdata, bool& insta);

void k2Search( FAM& FamG, vector<FAM>& DAG, vector<NODE>& Tree, int g, char op);

void searchDAGs( vector<FAM>& DAG, vector<NODE>& Tree, char op );

#define SEARCHDAGS_H_





#endif /* SEARCHDAGS_H_ */
