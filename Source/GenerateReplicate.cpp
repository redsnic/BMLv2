/*
GenerateReplicate.cpp
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
*/

#include "GenerateReplicate.h"

/*	Reads the input data file, identifies taxon & gene names and initializes the taxon data in variable "Tree".					*/

/**
 * Prepares the DAG used for output
 *
 * @param DAG             input   (a bayesan network)
 * @param SimDAG          output  (a DAG with probability parameters)
 * @param GeneLabels      name of genes
 * @param char_label      labels associated to nodes
 */
void parseDAG(vector<FAM>& DAG, vector<SimVar>& SimDAG, vector<string>& GeneLabels, vector< vector<int> >& char_label)
{
    SimDAG.resize(DAG.size());
    for(unsigned int i=0;i<DAG.size();i++)             // Initialization for each gene
    {
        SimVar temp;

        //temp.name=GeneLabels[char_label[i][0]];        //TODO print complete labels

        temp.name = "(";

        for(int c=0; c<char_label[i].size() ;c++){
        	temp.name += " ";
        	if(c == char_label[i].size()-1)
        		temp.name += GeneLabels[char_label[i][c]];
        	else
        		temp.name += GeneLabels[char_label[i][c]] + ",";
        }

        temp.name += " )";

        int ord=DAG[i].order;
        int numberOfParents=DAG[i].numberOfParents;
        for(unsigned int j=0;j<DAG.size();j++)
        {
            if(DAG[i].parents[j]==true){
                temp.parents.push_back(DAG[j].order);  // take parents informations
            }
        }
        int nparam=int(pow(2.0,numberOfParents+1));
        for(int j=0;j<nparam;j++)                      // and probability parameters
        {
            temp.k2Parameters.push_back(DAG[i].k2Parameters[j]);
        }
        temp.geneIndex=i;
        SimDAG[ord]=temp;
    }
}


/**
 * First stage in the phylogenetic tree preparation: vertex are associated with relative mutations.
 * "tree" will become a phylogenetic tree containig for each node (sample) the list of mutation for each gene
 * (so basically a row of the input matrix).
 *
 * @param tree      used for output, content is replaced
 * @param data      the input matrix
 */
void generateReplicate(vector<NODE>& tree, vector<vector<bool> >& data) {

	NODE vertex;
	vertex.parent = 0;
	vertex.numberOfDescendantSamples = 1;

	for (unsigned int i = 0; i < data[0].size(); i++) {    // the root is the normal genotype and has no mutations
		vertex.genotype.push_back(false);
	}

	tree.push_back(vertex);                                // Root initialized

	for (unsigned int i = 0; i < data.size(); i++) {       // the other vertex are initialized from the matrix
		vertex.genotype = data[i];
		tree.push_back(vertex);
	}
	tree[0].numberOfDescendantSamples = data.size();       // every node has the root as ancestor
}


/**
 * This procedure randomly creates nodes of an hypothetical output DAG
 *
 * @param tree           phylogenetic tree, will be randomly created (but maintaining its logical sense)
 * @param SimDAG         the output DAG
 * @param inputData      the data matrix given in input
 */
void  simulateDAG( vector<NODE>& tree, vector<SimVar>& SimDAG, vector< vector<bool> >& inputData)
{
    
	unsigned int sampleIndex=0;

    NODE v;
    v.parent=0;
    v.numberOfDescendantSamples=0;
    
	for(unsigned int i=0;i<SimDAG.size();i++)
    {
        v.genotype.push_back(false);
	}

	tree.push_back(v);

    while(sampleIndex<inputData.size())
	{

		vector<bool> simulatedGenotype(SimDAG.size());
		
        int numberOfHypotheticalMutations=0;
        
		for(unsigned int j=0;j<SimDAG.size();j++)
		{
			double x;
            do{
                x=(double) rand()/RAND_MAX;
            }while(x>=1);
            int y1=0;
            int numberOfParents=SimDAG[j].parents.size();
            double z=0.0;
            for(int k=0;k<numberOfParents;k++)
            {
                y1+=simulatedGenotype[SimDAG[j].parents[k]]*(int(pow(2.0,k)));        // parents state encoding
            }
            z=SimDAG[j].k2Parameters[y1]/(SimDAG[j].k2Parameters[y1]+SimDAG[j].k2Parameters[y1+(int(pow(2.0,numberOfParents)))] );
            if(SimDAG[j].k2Parameters[y1]+SimDAG[j].k2Parameters[y1+(int(pow(2.0,numberOfParents)))] <.999){cout<<"Prob norm \n";}
            if(SimDAG[j].k2Parameters[y1]+SimDAG[j].k2Parameters[y1+(int(pow(2.0,numberOfParents)))] >1.001){cout<<"Prob norm 2 \n";}
            if(z>x)
            {
                simulatedGenotype[j]=false;
            }else
            {
                simulatedGenotype[j]=true;
                numberOfHypotheticalMutations++;
            }
        }
        bool canAddNode=false;
        unsigned int sdt=0;

        while(canAddNode==false && sdt<inputData.size())
        {

            canAddNode=true;
            int numberOfMutations=0;
            for(unsigned int gch=0;gch<SimDAG.size();gch++)
            {
                if(simulatedGenotype[gch]==false && inputData[sdt][SimDAG[gch].geneIndex]==true){canAddNode=false;} // a mutation can't disappear
                if(inputData[sdt][SimDAG[gch].geneIndex]==true){numberOfMutations++;}
            }
            sdt++;
            if(numberOfMutations < 1){canAddNode=false;}                                                            // there must be mutations
        }

        if(canAddNode ==true && numberOfHypotheticalMutations >= 1)
        {
            sampleIndex++;
            for(unsigned int j=0;j<SimDAG.size();j++)
            {
                v.genotype=simulatedGenotype;
            }
            tree.push_back(v);                          // add node to the phylogenetic tree
        }
    }
}
