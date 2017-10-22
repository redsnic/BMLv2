/*
 *  ParseDataMatrix.cpp
 *  
 *
 Copyright (c) 2012-2014 Navodit Misra.
 
 BML is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 BML is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BML. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */

/*	Reads the input data file, identifies taxon & gene names and initializes the taxon data in variable "Tree".					*/

#include "ParseDataMatrix.h"

using namespace std; 

/**
 * Scans the data matrix to obtain the cutoff (over the number of mutations
 * of a single gene in the diffenrent samples) value necessary to have a number of genes
 * lower then the the number of samples (as it was done in the BML pubblication).
 *
 * @param inputMatrix    input matrix (without cutoff)
 * @param genesLabels    genes' Hugo codes
 */
void prepareCutoff(vector<vector<bool> >& inputMatrix, vector<string>& genesLabels, int forcedCutoff, bool autoCutoff){
	vector<int> counter(inputMatrix[0].size());

	for(unsigned int j = 0; j<inputMatrix[0].size(); j++){        // count number of mutations for each gene
		for(unsigned int i = 0; i<inputMatrix.size(); i++){
			counter[j] += inputMatrix[i][j];
		}
	}

	vector<int> numberOfValidGenesWithGivenCutoff(inputMatrix.size()+1);

	for(unsigned int i = 0; i<counter.size() ; i++){                        // count number of genes with a certain number of mutations
		numberOfValidGenesWithGivenCutoff[counter[i]]++;
	}

	for(int i = numberOfValidGenesWithGivenCutoff.size()-2; i>=0 ; i--){    // count number of genes with at most a certain number of mutations
		numberOfValidGenesWithGivenCutoff[i]+=numberOfValidGenesWithGivenCutoff[i+1];
	}

	int cutoff = 0;

	for(unsigned int i = 0; i<counter.size(); i++){                         // get cutoff value
		if(numberOfValidGenesWithGivenCutoff[i] < inputMatrix.size()){
			cutoff = i;
			break;
		}
	}

    if(!autoCutoff){         // use cutoff set by the user
		cutoff = forcedCutoff;
	}

	cout << "cutoff set with value " << cutoff << endl;

	transposeMatrix(inputMatrix);

	vector<vector<bool>> aux(inputMatrix.size());                            // init auxiliary data stractures
	for(unsigned int i = 0; i<inputMatrix.size(); i++){
		aux[i] = vector<bool>(inputMatrix[0].size());
	}

	vector<string> auxLabels(genesLabels.size());
	for(unsigned int i = 0; i < genesLabels.size(); i++ ){
		auxLabels[i] = "";
	}

	int validCount=0;

	for(unsigned int i = 0 ; i<inputMatrix.size(); i++){                     // update data matrix and
		if(counter[i]>=cutoff){                                               // genes' labels
			inputMatrix[i].swap(aux[validCount]);
			genesLabels[i].swap(auxLabels[validCount]);
			validCount++;
		}
	}

	inputMatrix.swap(aux);
	inputMatrix.resize(validCount);

	genesLabels.swap(auxLabels);                                              // finalize
	genesLabels.resize(validCount);

	transposeMatrix(inputMatrix);

}

/**
 * Standard counting sort to implement radixSort
 *
 * @param Data          the input boolean matrix
 * @param geneLabels    the genes' names corresponding to each column
 * @param row           the considered row for sorting
 */
void countingSort(vector< vector<bool> >& Data, vector<string>& geneLabels, int col){

	vector< vector<bool>> aux(Data.size());

	for(unsigned int i=0; i<aux.size(); i++){        // initializing auxiliary matrix
		aux[i] = vector<bool>(Data[0].size());
	}

	vector< string > aux_labels(geneLabels.size());
	int count[2] = {0,0};                            // counting

	for(unsigned int i = 0; i<Data.size() ; i++){    // counting for the given column
		count[Data[i][col]]++;
	}

	count[1] += count[0];                            // allocation vector

	for(int i = Data.size()-1; i>=0 ;i--){

		int elem = Data[i][col];
		int index = count[elem]-1;

		aux[index].swap(Data[i]);                    // swapping row

		aux_labels[index].swap(geneLabels[i]);       // swapping names
		count[elem]--;
	}

	Data.swap(aux);
	geneLabels.swap(aux_labels);

}

/**
 * modifies a boolean matrix to its transpose
 *
 * @param matrix (must be initialized)
 */
void transposeMatrix(vector< vector<bool> >& matrix){

	bool verticalGreaterSize = false;
	bool equalSizes = false;

	unsigned int oldSize=min(matrix.size(), matrix[0].size());

	cout << matrix.size() << " " << matrix[0].size() << endl;

	if(matrix.size() > matrix[0].size()){              // resize matrix
		for(unsigned int i = 0; i<matrix.size(); i++){
			matrix[i].resize(matrix.size());
		}
		verticalGreaterSize = true;
	}else if(matrix.size() < matrix[0].size()){
		for(unsigned int i = matrix.size(); i<matrix[0].size(); i++){
			matrix.push_back(vector<bool>(matrix[0].size()));
		}
	}else{
		equalSizes=true;
	}


	for(unsigned int i=0; i<matrix.size(); i++){       // transpose
		for(unsigned int j=i+1; j<matrix[i].size(); j++){
			bool temp = matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = temp;
		}
	}


	if( verticalGreaterSize && !equalSizes){           // resize matrix again
		matrix.resize(oldSize);
	}else if(!equalSizes){
		for(unsigned int i = 0; i<matrix.size(); i++){
			matrix[i].resize(oldSize);
		}
	}

	cout << matrix.size() << " " << matrix[0].size() << endl;

}

/**
 * Sorts each column of the matrix like it was a binary number
 * updating the order in the gene names to be consistent.
 *
 * @param Data          the input boolean matrix
 * @param geneLabels    the genes' names corresponding to each column
 */
void radixSort(vector< vector<bool> >& Data, vector<string>& geneLabels){

	transposeMatrix(Data);

	for(unsigned int i = 0; i<Data[0].size(); i++){
		countingSort(Data, geneLabels, i);
	}

	transposeMatrix(Data);
}

/**
 * Checks if two columns of a boolean matrix are equal
 *
 * @param matrix       boolean matrix
 * @param col1
 * @param col2
 * @return             true if the columns are equal, false otherwise
 */
bool equalCols(vector< vector<bool>>& matrix, int col1, int col2){
	for(unsigned int i=0; i<matrix.size(); i++){
		if(matrix[i][col1] != matrix[i][col2])
			return false;
	}
	return true;
}

/**
 * Associates together the columns with exactly the same mutations in all the samples
 * @param Data          rows of the input matrix
 * @param char_label    used for output, content is replaced
 * @return false if the input matrix is of all zeroes
 */
bool eliminateRedundantCharacters(vector< vector<bool> >& Data, vector<vector<int> >& char_label, vector<string>& geneLabels)//discard repeated columns
{	

	cout<< "Grouping columns with exactly the same mutations in all the samples:" << endl;

    if(Data.size()==0 || Data[0].size()==0){return false;}

	/* initialization of shouldThisColumnBeErased (like a row) */

	radixSort(Data, geneLabels);

	int originalSize = Data[0].size();

	vector< vector<bool>> aux(Data.size());

	for(unsigned int i=0; i<aux.size(); i++){
		aux[i] = vector<bool>(Data[i].size());
	}

	int numberOfUniqueCols = 0;

	for(unsigned int col = 0; col < Data[0].size(); col++){
		vector<int> groupReference;
		groupReference.push_back(col);

		for(unsigned int row = 0; row<Data.size(); row++){
			aux[row][numberOfUniqueCols] = Data[row][col];
		}
		numberOfUniqueCols++;

		int reference = col;
		col++;
		while(col < Data[0].size() && equalCols(Data, reference, col)){
			groupReference.push_back(col);
			col++;
		}
		char_label.push_back(groupReference);
		col--;
	}

	for(unsigned int i=0; i<aux.size(); i++){
		aux[i].resize(numberOfUniqueCols);
	}

	Data.swap(aux);

	/* debug output */

	cout << endl;
	cout << endl;

	for(unsigned int i = 0; i<Data.size(); i++){
		for(unsigned int j = 0; j<Data[0].size(); j++){
			if(Data[i][j])
				cout<< "X" <<" ";
			else
				cout << "  ";
		}
		cout << endl;
	}


	//TODO memory

	cout<< "Number of columns reduced from " << originalSize <<" to "<<char_label.size()<<'\n';

    return true;
}

/**
 * Parse a string (inp) to an integer (num) if possible
 * @param inp input string
 * @param num output value (modifies its value) NOTE: if inp is not a number the value of num has no meaning
 * @return true if inp is a string containing an integer
 */
bool isInteger(string& inp, int& num) {
	num = 0;
	for (unsigned int i = 0; i < inp.length(); i++) {
		if ((int(inp[i]) - int('0') >= 0) && (int(inp[i]) - int('0') <= 9)) {
			num = num * 10 + int(inp[i]) - int('0');
		} else {
			return false;
		}
	}
	return true;
}


/**
 * Reads the data file, stores taxon and gene names as vectors and partially initializes the variable "Tree".
 * @param Data              used for output, content is replaced (will store the rows of the input matrix, so that the i-th element
 *                                                                will be true iff there is a mutation on the i-th gene)
 * @param TaxonLabels		used for output, content is replaced (will store the Taxa's names)
 * @param GeneLabels		used for output, content is replaced (will store the Genes' names)
 * @param inputFile			data matrix file to be read
 */
void parseDataMatrix(vector<vector<bool> >& Data, vector<string>& TaxonLabels, vector<string>& GeneLabels, string& inputFile)
{
	ifstream DataFile(inputFile.c_str());

	/* opening input file */
	if(!DataFile){throw 0;}else{cout<<"Reading "<<inputFile<<'\n';}


	/* reading the first two elements (with error handling */
	int NTaxa, NGenes;
	string inp;


	DataFile>>inp;
	if(!isInteger(inp, NTaxa))
	{
		cerr<<"\nNon integral entry "<< inp <<" for # taxa in file "<< inputFile <<'\n';
		throw 1;
	}
	DataFile>>inp;
	if(!isInteger(inp, NGenes))
	{
		cerr<<"\nNon integral entry "<< inp <<" for # genes/characters in file "<< inputFile <<'\n';
		throw 1;
	}

	/* storing gene names in an array for further use */

	for(int i=0;i<NGenes;i++)
	{
		DataFile>>inp;
		GeneLabels.push_back(inp);
	}

	/* hist stores an histogram dividing taxa as a function of the number of observed mutations */

	vector<int> hist(HISTOGRAM_SIZE);
    for(int i=0;i<HISTOGRAM_SIZE ;i++)   //  initialization of hist
    {
        hist[i]=0;
    }

    /* storing taxon or sample names in an array for further use (with input controls) */

	for (int i = 0; i < NTaxa; i++) {
		DataFile >> inp;
		TaxonLabels.push_back(inp);
		vector<bool> v;
		int numberOfMutations = 0;
		for (int j = 0; j < NGenes; j++) {
			string state;
			DataFile >> state;
			if (state == "0") {
				v.push_back(false);
			} else {
				if (state == "1") {
					v.push_back(true);
					numberOfMutations++;
				} else {
					cerr
							<< "\nProblem parsing the data matrix entry for taxon (row) # "
							<< i << " and gene (column) # " << j << '\n'
							<< "Either the entries are not 0/1 or Num Genes is incorrectly specified in the header as "
							<< NGenes << '\n';
					throw 1;
				}
			}
		}
		if (numberOfMutations < HISTOGRAM_SIZE) {        
			hist[numberOfMutations]++;
		} else {
			hist[HISTOGRAM_SIZE-1]++;
		}
		Data.push_back(v);
	}

	/* show histogram */

	cout << endl << "This table shows the number of taxa (second column) with a certain amount of mutated genes (first column)" << endl << endl;


    for(int i=0; i<HISTOGRAM_SIZE; i++)
    {
    	if(i == HISTOGRAM_SIZE -1){
        	cout << "> " <<i<<'\t'<<hist[i]<<endl;
        } else {
        	cout<<i<<'\t'<<hist[i]<<endl;
        }

    }

	cout<<"\nTaxa initialized\n"<<Data.size()<<" taxa and "<<NGenes<<" genes/characters found."<<'\n';
	DataFile.close();

}

