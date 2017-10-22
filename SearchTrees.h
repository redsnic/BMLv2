/*
 * SearchTrees.h
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

#ifndef SEARCHTREES_H_

#ifndef DS
#include "DataStructures.h"
#define DS
#endif

#include "Utils.h"
#include "SearchDAGs.h"
#include "InitTree.h"
#include "GenerateReplicate.h"


void relabelOutlier(vector<NODE>& Tree,  int gn, int pdad, bool& tnew);
void swapnodes(vector<NODE>& Tree, vector<FAM>& DAG, int& i, int& gramps, int& spart, bool& change);
void searchTrees( vector<FAM>& DAG, vector<NODE>& Tree);

#define SEARCHTREES_H_

#endif /* SEARCHTREES_H_ */
