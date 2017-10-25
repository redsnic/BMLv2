/*
 *  StructLearn.cpp
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

#include "StructLearn.h"
using namespace std; 

mutex mtx;

/**
 * Parallel execution of structLearn (refer to that procedure)
 * NOTE: oldScore and BayesanNetwork are shared and modified so any access to them must be protected
 */
void structLearnThreadExecution( vector<FAM>& bayesanNetwork, double& oldScore, int& numberOfTreesToGenerate, vector<vector<bool>>& data ,string& job){

    for(int seed=0; seed<numberOfTreesToGenerate; seed++)
    {
        /*	Define the basic network object		*/

        vector<FAM> dag(0);
        vector<NODE> phy(0);

        generateReplicate(phy,data);
        globalPruning(dag,phy,'n',job);
        /*	Initialize the missing data values in the                          *
         *	variable "Phy", cf. InitTree.cpp and perform local pruning		   */
        
        initTree(phy,dag);

        cout << "tree initialized\n seed " << seed+1 << " out of " << numberOfTreesToGenerate << "\n";

        localPruning(dag,phy);
        
        /* score computation */

        double newScore=0.0;
        
        /*	Perform local search operations in tree space for imputing missing data values		*
         *	in variable "tree" using the MPE (Most Probable Explanation) criterion. 			*/
        
        searchTrees(dag,phy);

        for(unsigned int i=0;i<dag.size();i++)
        {
            newScore+=dag[i].score-dag[i].numberOfMutationsProbability;    // compute score
        }

        cout << "New score: "<< newScore << endl;

        mtx.lock();

        // begin of the critical section

        if(newScore>oldScore)
        {
            cout << "NEW BEST SCORE: " << newScore << "(old: "<<oldScore<<")"<< endl;
            for(unsigned int i=0;i<bayesanNetwork.size();i++)              // update the bayesan network (= should be enough)
            {
                bayesanNetwork[i]=dag[i];
            }
            oldScore=newScore;                                             // had no sense putting this in the loop
        }
        
        // end of the critical section

        mtx.unlock();

    }
}


/**
 * Learns the bayesan network from the input data
 *
 * @param inputMatrix   the input data matrix
 * @param simDAG        used for output, content is replaced
 * @param inputFile     the input matrix
 * @param nTrees        number of random restarts
 * @param pthres  		threshold for inferring paths
 * @param job           the name given to the job
 * @return false if the data matrix is of zeroes only
 */
bool structLearn(vector< vector<bool> >& inputMatrix, vector<SimVar>& simDAG, string& inputFile, int& nTrees, double& pthres, string& job, int& numberOfThreads, int cutoff, bool autoCutoff)
{

    vector<string> taxonLabels;          // matrix data containers
    vector<string> geneLabels;
    vector<vector<int> > char_label;

    /* Read input file */

    parseDataMatrix(inputMatrix, taxonLabels, geneLabels, inputFile);

    prepareCutoff(inputMatrix, geneLabels, cutoff, autoCutoff);

    if (eliminateRedundantCharacters(inputMatrix, char_label, geneLabels) == false){
        return false;
    }

    vector<NODE> phylogeneticTree;
    vector<FAM> bayesanNetwork;

    generateReplicate(phylogeneticTree, inputMatrix);

    globalPruning(bayesanNetwork,phylogeneticTree,'y',job);

    /*	Initiate a subroutine for randomly reseeding the search procedure. Each tree serves as	*
     *	a starting point for a local search heuristic through tree space						*/

    double oldScore=INT_MIN;

    /* multi-threaded execution (reeseeding trees and maximizing scores) */

    int numberOfRepetitions = nTrees/numberOfThreads;
    int numberOfRepetitionsPlusOne = numberOfRepetitions+1;
    int extraRepetitions = nTrees%numberOfThreads;

    vector<thread> threads;

    for(int i=0; i<numberOfThreads; i++){
    	if(extraRepetitions>0){
    		threads.push_back(
                    thread( ref(structLearnThreadExecution) , ref(bayesanNetwork), ref(oldScore), ref(numberOfRepetitionsPlusOne), ref(inputMatrix), ref(job) ));
    		extraRepetitions--;
    	}else{
    		threads.push_back(
                    thread( ref(structLearnThreadExecution) , ref(bayesanNetwork), ref(oldScore), ref(numberOfRepetitions), ref(inputMatrix), ref(job) ));
    	}
    }

    for(int i=0; i<numberOfThreads; i++){
    	threads[i].join();
    }

    cout << "FINAL SCORE:" << oldScore << endl;

    /* preparing for output */

    parseDAG(bayesanNetwork, simDAG, geneLabels,char_label);
    writeLandscape(simDAG,phylogeneticTree, pthres, job);
    /*                      Parametric Bootstrap                                                */
    


    return true;
}

/**
 * Parallel execution of bootstrap replicates analysis (refer to bootstrapAnalysis)
 */
void bootstrapAnalysisThreadExecution(vector< vector<bool> >& Data, vector<SimVar>& SimDAG, int& NRep, string& job,
		vector<vector<double> >& landrel1, vector<vector<double> >& landrel2, vector< vector<int> >& netrel1, vector< vector<int> >& netrel2,
		vector< vector<string> >& edgeList, vector< vector<double> >& edgeProbability){


    for(int BootstrapReplication=1; BootstrapReplication<NRep+1;BootstrapReplication++)  // parallel Execution
    {
        vector<NODE> tree;
        vector<FAM> bayesNP;
        
        cout<<"starting simulation "<<BootstrapReplication<<" out of "<<NRep<<'\n';

        simulateDAG(tree, SimDAG, Data);         // compute new tree
        
        cout<<"\nGenerated bootstrap replicate "<<BootstrapReplication<<" out of "<<NRep<<'\n';

        globalPruning(bayesNP,tree,'n',job);

        searchDAGs(bayesNP,tree,'n');

        for(unsigned int i=0;i<netrel2.size();i++)
        {
            for(unsigned int j=0;j<netrel2.size();j++)
            {
                netrel2[i][j]+=int(bayesNP[i].parents[j])+int(bayesNP[j].parents[i]);    //compute netrel2
            }
        }


        landBoot(SimDAG, landrel1, bayesNP);           // compute probability of Tree+OBSDAG edges
        
        /*	Parse each input tree and assign the missing data values in the						*
         *	variable "Tree", cf. InitTree.cpp.													*/
        
        vector<FAM> DAG(bayesNP.size());					/*	Define the basic network object	*/
        
        initTree(tree, bayesNP);

        initDAG(DAG, bayesNP);

        localPruning(DAG,tree);


        /*	Perform local search operations in tree space for imputing missing data values		*
         *	in variable "Tree" using the MPE (Most Probable Explanation) criterion. 			*/
        
        searchTrees(DAG,tree);


        cout<< BootstrapReplication << " out of " << NRep << " bootstrap replicates analyzed\n";
        for(unsigned int i=0;i<netrel1.size();i++)
        {
            for(unsigned int j=0;j<netrel1.size();j++)
            {
                netrel1[i][j]+=int(DAG[i].parents[j])+int(DAG[j].parents[i]);  // 0-1-2 values
            }
        }

        landBoot(SimDAG, landrel2, DAG);              // compute probability of Tree+OBSDAG edges

        edgeBoot(edgeList,SimDAG,edgeProbability, DAG);
    }


}


/**
 * This procedure generates the bootstrap replicates and repeats
 * the analysis to compute the confidence intervals of probabilities
 *
 * @param Data             input data matrix
 * @param SimDAG
 * @param NRep             number of bootrstrap replications
 * @param job              name of the job
 * @param numberOfThreads
 */
void bootstrapAnalysis(vector< vector<bool> >& Data, vector<SimVar>& SimDAG, int& NRep, string& job, int numberOfThreads)
{

    /* Edges present in Tree+OBSDAG, OBSDAG and SimDAG resp. */
    vector< vector<int> > netrel1, netrel2, netrel3;
    vector< vector<string> > edgeList;
    initNetRel(netrel1, netrel2, netrel3, edgeList, SimDAG);

    edgeList.erase(edgeList.begin(), edgeList.end());

    /* Bootstrap probabilities for the top ordered genes for Tree+OBSDAG and OBSDAG (not to be confused with the ordering of mutations) */
    vector<vector<double> > landrel1(SimDAG.size());
    vector<vector<double> > landrel2(SimDAG.size());

    /* Bootstrap probabilities for each pair of genes sharing an edge, 2 single mutant states and the double mutant genotype using Tree+OBSDAG */
    vector< vector<double> > edgeProbability;

    /* prepare threads vector */

    int numberOfExecutionsPerThread = NRep/numberOfThreads;   // number of executions for each thread
    int extraExecutions = NRep%numberOfThreads;               // numberOfExecutionsPerThread can not be exactly NRep
    cout << "Level of multi-threading: " << numberOfThreads <<"\n";
    vector<thread> threads;

    /* vectors needed for multi-threading */

    vector<vector< vector<int> >> localNetrels1, localNetrels2, localNetrels3;
    vector<vector< vector<string> >> localEdgeLists;
    vector<vector<vector<double> >> localLandrels1, localLandrels2;
	vector<vector< vector<double>>> localEdgeProbabilities;

    for(int i=0; i<numberOfThreads; i++ ){

    	vector< vector<int> > tempNetrel1, tempNetrel2, tempNetrel3;
    	vector< vector<string> > tempEdgeList;

    	vector<vector<double> > tempLocalLandrels1 (SimDAG.size());
		vector<vector<double> > tempLocalLandrels2 (SimDAG.size());

    	initNetRel(tempNetrel1, tempNetrel2, tempNetrel3, tempEdgeList, SimDAG);
    	localNetrels1.push_back(tempNetrel1);
    	localNetrels2.push_back(tempNetrel2);
    	localNetrels3.push_back(tempNetrel3);
    	localEdgeLists.push_back(tempEdgeList);

    	vector< vector<double> > tempEdgeProbabilities (3*tempEdgeList.size());

    	localLandrels1.push_back(tempLocalLandrels1);
    	localLandrels2.push_back(tempLocalLandrels2);
    	localEdgeProbabilities.push_back(tempEdgeProbabilities);

    }

    /* parallel execution */

    int numberOfExecutionsPerThreadPlusOne = numberOfExecutionsPerThread+1;

    for(int i=0; i<numberOfThreads; i++ ){
    	cout << "launching thread: " << i <<endl;
    	if(extraExecutions>0){

			threads.push_back(
					thread(ref(bootstrapAnalysisThreadExecution), ref(Data), ref(SimDAG),
							ref(numberOfExecutionsPerThreadPlusOne), ref(job), ref(localLandrels1[i]),
							ref(localLandrels2[i]), ref(localNetrels1[i]), ref(localNetrels2[i]), ref(localEdgeLists[i]),
							ref(localEdgeProbabilities[i])));
    		extraExecutions--;
    	} else {
			threads.push_back(
					thread(ref(bootstrapAnalysisThreadExecution), ref(Data), ref(SimDAG),
							ref(numberOfExecutionsPerThread), ref(job), ref(localLandrels1[i]),
							ref(localLandrels2[i]), ref(localNetrels1[i]), ref(localNetrels2[i]), ref(localEdgeLists[i]),
							ref(localEdgeProbabilities[i])));
    	}
    }

    for(int i=0; i<numberOfThreads; i++ ){
    	threads[i].join();
    }

    cout << "merging data from threads..." << endl;

    /* merge data */

    for(int thr=0; thr<min(numberOfThreads, NRep); thr++ ){

    	//cout << netrel1.size() << " " << netrel2.size() << " " << netrel3.size() << endl;


        for(unsigned int i=0; i<netrel1.size() ;i++)
        {
            for(unsigned int j=0;j<netrel1.size();j++)
            {
                netrel1[i][j]+= localNetrels1[thr][i][j];
                netrel2[i][j]+= localNetrels2[thr][i][j];
            }
        }

        for(unsigned int i = 0; i<localLandrels1[thr].size(); i++){
        	landrel1[i].insert(landrel1[i].end(), localLandrels1[thr][i].begin() , localLandrels1[thr][i].end());
        	landrel2[i].insert(landrel2[i].end(), localLandrels2[thr][i].begin() , localLandrels2[thr][i].end());
        }


        edgeProbability.insert(edgeProbability.end(), localEdgeProbabilities[thr].begin() , localEdgeProbabilities[thr].end());

        for(unsigned int i = 0; i<localEdgeLists[thr].size(); i++){
        	edgeList.push_back(localEdgeLists[thr][i]);
        }

    }

    cout << "finalizing ...";

    /* final stage */

    tpfpAnalysis(netrel1,netrel2,netrel3,SimDAG,NRep,job);

    cout << "...";

    probBootstrapAnalysis(landrel1,landrel2,SimDAG,job); // TODO correct names in sub routine

    cout << "...";

    probEdBoot(edgeList, edgeProbability,job);

    cout << " done!" << endl;
}
