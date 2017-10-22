/*
 * StructLearn.h
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
#include <thread>
#include <mutex>

#ifndef STRUCTLEARN_H_

#ifndef DS
#include "DataStructures.h"
#define DS
#endif

#include "IOHelper.h"
#include "ParseDataMatrix.h"
#include "Utils.h"
#include "InitTree.h"
#include "SearchTrees.h"
#include "mainmenu.h"
#include "ui_mainmenu.h"

bool structLearn(vector< vector<bool> >& Data, vector<SimVar>& SimDAG, string& inputFile, int& NTrees, double& pthres, string& job, int& numberOfThreads, int cutoff, bool autoCutoff);

void bootstrapAnalysis(vector< vector<bool> >& Data, vector<SimVar>& SimDAG, int& NRep, string& job, int numberOfThreads);



#define STRUCTLEARN_H_


#endif /* STRUCTLEARN_H_ */
