/*
 * GenerateReplicate.h
 *
 *  Created on: 24 ago 2017
 *      Author: rossi
 */

#ifndef GENERATEREPLICATE_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>

#ifndef DS
#include "DataStructures.h"
#define DS
#endif
#include "Utils.h"
#include "SearchDAGs.h"

void parseDAG(vector<FAM>& DAG, vector<SimVar>& SimDAG, vector<string>& GeneLabels, vector< vector<int> >& char_label);
void generateReplicate(vector<NODE>& Tree, vector<vector<bool> >& Data);
void  simulateDAG( vector<NODE>& Tree, vector<SimVar>& SimDAG, vector< vector<bool> >& sdata);
#define GENERATEREPLICATE_H_

#endif /* GENERATEREPLICATE_H_ */
