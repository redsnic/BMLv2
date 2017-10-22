/*
 * MAFtoMatrixTranslator.h
 *
 * utility to get informations from MAF (multiple annotation format) files
 * and obtain the data matrix used by BML.
 *
 *  Created on: 03 set 2017
 *      Author: Nicolo' Rossi
 */

#ifndef MAFTOMATRIXTRANSLATOR_H_

#include <iostream>
#include <string>
#include <fstream>
#include<vector>
#include <string.h>
#include <algorithm>
#include <StringTokenizer.h>

#define MAFTOMATRIXTRANSLATOR_H_

using namespace std;

namespace std {


class MAFtoMatrixTranslator {
public:
	MAFtoMatrixTranslator();
    bool read(string path, bool all, bool standard, vector<string>& mutations, bool useCustomGenes, vector<string>& genes);
	bool save(string path);
	virtual ~MAFtoMatrixTranslator();

private:
	vector<string> genesList;
	vector<string> samplesList;
	vector<vector<bool> > dataMatrix;   // [sample][gene]
	int numberOfGenes;
	int numberOfSamples;
	bool prepareDataMatrix(ifstream& mafInput);
    bool populateDataMatrix(ifstream& mafInput, vector<string>& selectedTypesOfMutation, bool all, bool standard, bool useCustomGenes, vector<string>& genes);
	void skipComments(ifstream& mafInput);
	void sortAndEliminateRepetitions(vector<string>& v);
	int dichotomic_search(vector<string>& v, string value);
    void removeUnusedSamples();
};

} /* namespace std */

#endif /* MAFTOMATRIXTRANSLATOR_H_ */
