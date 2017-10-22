/*
 * Utils.h
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

#ifndef UTILS_H_
#include "mainmenu.h"
#include "ui_mainmenu.h"

#ifndef DS
#include "DataStructures.h"
#define DS
#endif

void initializeDAG(vector<FAM>& DAG);

bool likelihoodTestLocal(int i, int j,  vector<NODE>& Tree );

void initNumberOfMutationProbability(  vector<FAM>& BN, vector<NODE>& Phylogeny) ;

void localPruning( vector<FAM>& DAG, vector<NODE>& Tree ) ;

void initDAG(vector<FAM>& DAG, vector<FAM>& BN);

bool likelihoodTestGlobal(double& olap, double& Atot, double& Btot, double& NSam);

void globalPruning( vector<FAM>& BayesN, vector<NODE>& Tree, char op, string& job);

#define UTILS_H_



#endif /* UTILS_H_ */
