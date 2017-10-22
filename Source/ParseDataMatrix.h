/*
 * ParseDataMatrix.h
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

#ifndef PARSEDATAMATRIX_H_
#include "Utils.h"
#ifndef DS
#include "DataStructures.h"
#define DS
#endif


void delete_Col(vector< vector<bool> >& Data, int col_no);
bool eliminateRedundantCharacters(vector< vector<bool> >& Data, vector<vector<int> >& char_label, vector<string>& geneLabels); // TODO roll
bool isInteger(string& inp, int& Num);
void parseDataMatrix(vector<vector<bool> >& Data, vector<string>& TaxonLabels, vector<string>& GeneLabels, string& dname);
void countingSort(vector< vector<bool> >& Data, vector<string>& geneLabels, int row);
void radixSort(vector< vector<bool> >& Data, vector<string>& geneLabels);
void prepareCutoff(vector<vector<bool> >& inputMatrix, vector<string>& genesLabels, int forcedCutoff, bool autoCutoff);
void transposeMatrix(vector< vector<bool> >& matrix);

#define PARSEDATAMATRIX_H_

#define HISTOGRAM_SIZE 21

#endif /* PARSEDATAMATRIX_H_ */
