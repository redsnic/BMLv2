/*
 *  DataStructures.h
 *  
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
 */

#include "ui_mainmenu.h"

/*	Defines the basic data types used elsewhere in the code.								*/

using namespace std;

/*	A Tree is a vector of these NODES, it is also our data matrix.							*/
struct NODE {
	int parent;		    				 //	Parent node on the tree.
	vector<int> children;
	vector<int> numberOfMutationsOnEdge;	     //	# of mutations on the incoming edge.
	vector<bool> genotype;			     //	Binary genotype (mutated -> true).
	int numberOfDescendantSamples;       //	Number of samples (taxa) that are descendants.
};

/*	Local structure of the Bayes Net for each family (parents + param for each gene) (Markovian Blanket)	*/

struct FAM {
	vector<bool> parents;			      // if an element is true its index refers to a node's parent, like in an adjacency matrix
	vector<bool> potentialParents;
	vector<double> k2Parameters;          // parameters computed by K2  procedure
	double numberOfMutations;             // number of mutation for the given gene
	double numberOfMutationsProbability;  // log probability of having this number of mutations (assuming a Bernoulli distribution)
	double score;
	int numberOfParents;
	int order;                            //  Gene order during OBS
	bool hasNoParent;
};

/* a node of the Bayesan Network */

struct SimVar {
	string name;
	vector<int> parents;
	vector<double> k2Parameters;
	int geneIndex;
};


/* input structure needed to start execute as a thread */

struct Input {
    string job;
    string inputFile;
    int numberOfThreads;
    int nTrees;
    double pthres;
    bool performBootstrap;
    int nRep;
    bool autoCutoff;
    int cutoff;
};
