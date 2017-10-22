/*
 * IOHelper.h
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
#include <limits>
#include "MAFtoMatrixTranslator.h"

#ifndef IOHELPER_H_


#ifndef DS
#include "DataStructures.h"
#define DS
#endif
#include "Utils.h"

void getMatrixFromMAF();

void initNetRel(vector<vector<int> >& netrel1, vector<vector<int> >& netrel2,
		vector<vector<int> >& netrel3, vector<vector<string> >& EdList,
		vector<SimVar>& SimDAG);

void tpfpAnalysis(vector<vector<int> >& NetRel1, vector<vector<int> >& NetRel2,
		vector<vector<int> >& NetRel3, vector<SimVar>& SimDAG, int& NRep,
		string& job);
void configParam(vector<SimVar>& SimDAG, vector<bool>& config, double& pratio);
void writeLandscape(vector<SimVar>& SimDAG, vector<NODE>& tree, double& pthres, string& job);
void landBoot(vector<SimVar>& SimDAG, vector<vector<double> >& landrel,
		vector<FAM>& DAG);

void edgeBoot(vector<vector<string> >& EdList, vector<SimVar>& SimDAG,
		vector<vector<double> >& edrel, vector<FAM>& DAG);

void probEdBoot(vector<vector<string> >& EdList, vector<vector<double> >& edrel,
		string& job);

void probBootstrapAnalysis(vector<vector<double> >& NetRel1,
		vector<vector<double> >& NetRel2, vector<SimVar>& SimDAG, string& job);

void writeDAG(vector<FAM>& DAG, vector<string>& GeneLabels,
        vector<vector<int> >& char_label, char op);

void exceptionHandler(int& ex);

#define IOHELPER_H_



#endif /* IOHELPER_H_ */
