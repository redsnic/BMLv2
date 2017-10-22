/*
 * MAFtoMatrixTranslator.cpp
 *
 *  Created on: 03 set 2017
 *      Author: rossi
 */

#include "MAFtoMatrixTranslator.h"

namespace std {

MAFtoMatrixTranslator::MAFtoMatrixTranslator() {

	numberOfGenes = 0;
	numberOfSamples = 0;

}


/**
 * This procedure takes in input a path to a MAF (multiple annotations format)
 * file and translates it internally to generate a
 * gene/patient boolean matrix that has 'true' value in a cell iff the relative patient
 * had a mutation on the relative gene (so a BML input file).
 *
 * @param path     path to MAF file (data is only read, file is not modified)
 * @param all                        true if all mutations should be used
 * @param standard                   true if only the BML compatible mutations should be used
 * @param mutations                  list of mutations choosen by the user
 * @param useCustomGenes             true if the genes to be used were choosen by the user
 * @param genes                      list of mutations choosen by the user (unused if usedCustomGenes is false)
 * @return         true if translation was successful, false otherwise
 */
bool MAFtoMatrixTranslator::read(string path, bool all, bool standard, vector<string> &mutations, bool useCustomGenes, vector<string> &genes){

	ifstream mafInput (path.data());     // open input file
	if (mafInput.bad() || !mafInput.is_open()){
		mafInput.close();
		return false;
	}

	skipComments(mafInput);

	string mafFieldHeader;               // skip mafField header
	getline(mafInput, mafFieldHeader);

	skipComments(mafInput);

	streampos firstLinePosition = mafInput.tellg();    // remember matrix starting point

	if( !prepareDataMatrix(mafInput) ){  // initialize basic matrix structure
		mafInput.close();
		return false;
	}

	mafInput.clear();
	mafInput.seekg(firstLinePosition);   // go back to the beginning of the matrix

    if( !populateDataMatrix(mafInput, mutations, all, standard, useCustomGenes, genes) ){ // set matrix content
		mafInput.close();
		return false;
	}

	mafInput.close();

	return true;                         // initialization complete

}


/**
 * Prepares the boolean matrix used to generate the input file reading the data from the input stream
 *
 * @param mafInput                   the input stream
 * @param selectedTypesOfMutation    list of mutations choosen by the user
 * @param all                        true if all mutations should be used
 * @param standard                   true if only the BML compatible mutations should be used
 * @param useCustomGenes             true if the genes to be used were choosen by the user
 * @param genes                      list of mutations choosen by the user (unused if usedCustomGenes is false)
 * @return            true if the operation was successful
 */
bool MAFtoMatrixTranslator::populateDataMatrix(ifstream& mafInput, vector<string>& selectedTypesOfMutation, bool all, bool standard , bool useCustomGenes, vector<string>& genes){

    string line;

    if(useCustomGenes){
        this->sortAndEliminateRepetitions(genes);
        genesList=genes;
        numberOfGenes=genes.size();
    }

	int sampleIndex = 0, geneIndex = 0;

	/* data matrix initialization */

	for(int i = 0; i<numberOfSamples ; i++){
		vector<bool> v;
 		dataMatrix.push_back(v);
		for(int j = 0; j<numberOfGenes ; j++){
			dataMatrix[i].push_back(false);
		}
	}

	/* populating the data matrix */

	string oldSample = "";
	string type = "";
	string genere = ""; // SNP-CNV-DEL-INS

	try {
		while ( !mafInput.eof() ) {                         // read all the MAF file

            bool ignore = false;

			getline(mafInput, line);                        // and for each line (mutation)

            StringTokenizer tok(line, "\t");

            string field = tok.next();
			for (int i = 1; i<17 && !mafInput.eof() ; i++) {

				if (i == 1) {                               // get gene and sample name

                    if(useCustomGenes && genesList.end() == find(genesList.begin(), genesList.end(), field)){//!binary_search(genesList.begin(), genesList.end(), field, [](std::string a, std::string b){ return b < a; })){   // skip not choosen genes
                        ignore =true;
                        break;
                    }

                    if(field == "Unknown"){                  // do not consider unknown genes
                        ignore =true;
                        break;
					}else{

                        geneIndex =  dichotomic_search(genesList, field);      // mutated gene
					}
				}else{
                    field = tok.next();
                    //cout << field << endl;
				}

				if(i == 9){           // mutation type
                    type = field;
				}

				if(i==10){            // SNP - CNV
                    genere = field;
				}

				if (i == 16) {        // sample

                    if(field != oldSample){
                        sampleIndex = dichotomic_search(samplesList, field);  // considered sample
					}
                    oldSample = field;
				}
			}
            if (!ignore){
                bool found = false;
                if( all || standard ){
                    found = true;
                }else{
                    found = any_of(begin(selectedTypesOfMutation), end(selectedTypesOfMutation), [&](string i){ return (type == i); });
                }

                if((genere != "SNP" && genere != "DEL"  && genere != "INS") && standard){  // keep only SNP or INDEL (STD)
                    //cout << "Non SNP mutation removed (" << genere <<")"<< endl;
                    found = false;
                }

                if( (type == "Silent" || type == "silent") && standard){  // remove silent mutations (STD)
                    //cout << "Silent mutation removed" << endl;
                    found = false;
                }

                if(found){                                                    // if it is a valid mutation type
                    cout << sampleIndex << " : " << geneIndex << endl;
                    dataMatrix[sampleIndex][geneIndex] = true;                // set mutation
                }
            }
		}
	} catch (const int e) {
		return false;
	}

    removeUnusedSamples();

	return true;

}

/**
 * @brief MAFtoMatrixTranslator::removeUnusedSamples
 *
 * removes the unused samples from the data matrix
 * ( so the rows of the matrix of zeroes
 */
void MAFtoMatrixTranslator::removeUnusedSamples(){
    for(unsigned int i = 0; i<dataMatrix.size(); i++){
        bool found = false;
        for(unsigned int j=0; j<dataMatrix[i].size(); j++){
            if(dataMatrix[i][j] == true){
                cout << i << " : " << j << endl;
                found=true; break;
            }
        }
        if(!found){
            cout << i << endl;
            dataMatrix.erase(dataMatrix.begin()+i);
            samplesList.erase(samplesList.begin()+i);
            numberOfSamples--;
            i--;
        }
    }
}

/**
 * binary search for a string vector without repetitions
 *
 * @param v            must have no repetitions and be sorted
 * @param lowerBound   range in which to search
 * @param upperBound
 * @param value        value to search ( must be present )
 * @return             the position of value
 */
int MAFtoMatrixTranslator::dichotomic_search(vector<string>& v, string value){

	int lowerBound=0;
	int upperBound=v.size()-1;

	int middle = (lowerBound + upperBound )/2;

	do{
        if(v[middle] != value){
            if(value > v[middle]){
                lowerBound = middle+1;
                middle = (lowerBound + upperBound )/2;
            }else{
                upperBound = middle-1;
                middle = (lowerBound + upperBound )/2;
            }
        }
	}while(v[middle] != value);

	return middle;

}


/**
 * reads the maf file to get genes names and patients IDs (1st and 14th columns of MAF file)
 *
 * @param mafInput   stream of a maf file pointing to the first line of the input matrix
 * @return  true if operation was successful
 */
bool MAFtoMatrixTranslator::prepareDataMatrix(ifstream& mafInput){

	string gene, sample, line;

	int count  = 1;

	try {
		while (!mafInput.eof()) {                           // read all the MAF file



			getline(mafInput, line);                        // and for each line (mutation)
			count++;
            StringTokenizer tok(line, "\t");
            string field = tok.next();

            for (int i = 1; i<17 && !mafInput.eof() ; i++) {

                if (i == 1) {                               // get gene and sample name
                    if(field == "Unknown"){                  // remove unknown genes
						break;
					}
                    genesList.push_back(field);
				}else{
                    field = tok.next();
				}
                if (i == 16) {	
                    samplesList.push_back(field);
				}
			}
		}

	} catch (const int e) {
		return false;
	}

	sortAndEliminateRepetitions(genesList);
	sortAndEliminateRepetitions(samplesList);

	numberOfGenes = genesList.size();
	numberOfSamples = samplesList.size();

	return true;

}

/**
 * Sorts in lexicographical order a vector v of strings and eliminates repeated words
 *
 * @param v vector of strings
 */
void MAFtoMatrixTranslator::sortAndEliminateRepetitions(vector<string>& v){
	sort(v.begin(), v.end());

	vector<string> aux;
	if(v.size() >0){
		aux.push_back(v[0]);
    }

	for(unsigned int i = 1 ; i < v.size() ; i++){
		if(!(v[i-1] == v[i])){
			aux.push_back(v[i]);
		}
	}

	v=aux;


}
/**
 * Skips comments (lines with '#' as first non space character)
 * and white lines from an input stream.
 *
 * @param inputFile     input stream
 */
void MAFtoMatrixTranslator::skipComments(ifstream& inputFile){

	streampos oldPosition = inputFile.tellg();    // remember old position to avoid to skip a non-comment line
	string line;

	while(!inputFile.eof()){                      // for each line
		getline(inputFile, line);
		size_t firstCharacter = line.find_first_not_of(" \t");
		if( firstCharacter == string::npos || line[firstCharacter] == '#' ){  // if it is blank or a comment
			 oldPosition = inputFile.tellg();                                 // skip it and update stream position
		} else {
			 inputFile.seekg(oldPosition);                                    // else go back and return
			 return;
		}
	}

}

/**
 * saves the output data matrix in a given location
 *
 * @param path     where to create the output file
 * @return         true if the operation was successful
 */
bool MAFtoMatrixTranslator::save(string path){

	ofstream outputFile (path.data());

	if(outputFile.bad()){
		return false;
	}

	outputFile << numberOfSamples;	                   // header
	outputFile << "\t";
	outputFile << numberOfGenes;
	outputFile << endl;


	for(unsigned int i=0; i<genesList.size(); i++){    // genes
		outputFile << genesList[i] << "\t";
	}
	outputFile << endl;

	for(unsigned int i=0; i<samplesList.size(); i++){  // samples with genome
		outputFile << samplesList[i] << "\t";
		for(unsigned int j=0; j<genesList.size(); j++){
			outputFile << dataMatrix[i][j] << "\t";
		}
		outputFile << "\n";
	}

	outputFile.close();

	return true;

}

MAFtoMatrixTranslator::~MAFtoMatrixTranslator() {

}

} /* namespace std */
