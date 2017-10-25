/*
 *  Utils.cpp
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

#include "Utils.h"

using namespace std; 

/* Utilities for pruning and DAG initialization */

/**
 * Basic DAG initialization (from globalPruning output)
 *
 * @param DAG        the bayesan network to be initialized
 */
void initializeDAG(vector<FAM>& DAG)
{
    FAM famtemp;
    
    for(unsigned int i=0;i<DAG.size();i++)
    {
        famtemp.parents.push_back(false);
    }
    famtemp.numberOfParents=0;
    famtemp.numberOfMutations=0.0;
    famtemp.score=0.0;

    for(unsigned int i=0;i<DAG.size();i++)
    {
        DAG[i].parents=famtemp.parents;
        DAG[i].numberOfParents=famtemp.numberOfParents;         // 0
        DAG[i].score=famtemp.score;                             // 0
        DAG[i].numberOfMutations=famtemp.numberOfMutations;     // 0
    }
}

/**
 * Computes the probability of having the exactly the amount of mutations considered for each gene
 * and the probability of having a mutation in a gene. The data is then stored in the bayesan network.
 *
 * @param bayesanNet     the bayesan network.
 * 						 Will be updated with numberOfMutations for each gene,
 * 						 probability of having such number of mutations considering a Bernoulli distribution (so randomly)
 * 						 and the probability of having a mutation in a gene randomly considering a sample.
 *
 * @param phylogeny      the phylogeny (to analyze the number of mutations)
 */
void initNumberOfMutationProbability(vector<FAM>& bayesanNet, vector<NODE>& phylogeny)
{
    for(unsigned int i=0;i<bayesanNet.size();i++)                                // for each gene
    {
        double numberOfMutations=0.0;
        double numberOfSamples=phylogeny.size();
        for(unsigned int k=0;k<phylogeny.size();k++)
        {
            if(phylogeny[k].genotype[i]==true){numberOfMutations+=1.0;}          // count the number of mutations
        }

        bayesanNet[i].numberOfMutations=numberOfMutations;                       // store it
        double theta=numberOfMutations/numberOfSamples;                          // compute its probability
        if(theta>0.0)
        {
            bayesanNet[i].numberOfMutationsProbability = - numberOfMutations*log(theta)-(numberOfSamples-numberOfMutations)*log(1.0-theta);
            // log probability of i ( Bernoulli distribution )
        }
        else
        {
            bayesanNet[i].numberOfMutationsProbability=0.0;
        }

        bayesanNet[i].k2Parameters.erase(bayesanNet[i].k2Parameters.begin(),bayesanNet[i].k2Parameters.end());
        bayesanNet[i].k2Parameters.push_back(1.0-theta);                          // k2Parameters contains now the probability of a mutation in the given gene
        bayesanNet[i].k2Parameters.push_back(theta);                              // and its complement (first parameter)
    }
}

/**
 * Computes likelihood score and compares it to logarithm of the total number of samples
 * (refer to chapter 2.2 (formula 9) of the supplementary informations)
 * (local pruning)
 * @param i        potential parent (gene)
 * @param j        potential child  (gene)
 * @param tree     phylogenetic tree
 * @return         true if the test is passed
 */
bool likelihoodTestLocal(int i, int j,  vector<NODE>& Tree )
{
    double cdat[2][2];
    cdat[0][0]=0;
    cdat[0][1]=0;
    cdat[1][0]=0;
    cdat[1][1]=0;
    double numberOfSamples=Tree.size();
    double MInf=0.0;

    double nmin2=0.0;
    for(unsigned int k=0;k<Tree.size();k++)  // for each sample
    {	
        cdat[Tree[k].genotype[i]][Tree[k].genotype[j]]++;  // count the number of times i and j were in a certain respective state
                                                           // (i, j both mutated, only j mutated, only i mutated and none of them mutated)
    }
    
    MInf=cdat[0][0]*log((numberOfSamples*cdat[0][0])/((cdat[0][0]+cdat[0][1])*(cdat[0][0]+cdat[1][0])));
    if(cdat[0][1]!=0)
    {
        MInf+= cdat[0][1]*log((numberOfSamples*cdat[0][1])/((cdat[0][0]+cdat[0][1])*(cdat[0][1]+cdat[1][1])));
    }
    if(cdat[1][0]!=0)
    {
        MInf+= cdat[1][0]*log((numberOfSamples*cdat[1][0])/((cdat[1][0]+cdat[1][1])*(cdat[0][0]+cdat[1][0])));
        
    }
    if(cdat[1][1]!=0)
    {
        MInf+=cdat[1][1]*log((numberOfSamples*cdat[1][1])/((cdat[1][0]+cdat[1][1])*(cdat[0][1]+cdat[1][1])));
        nmin2=cdat[1][1]/(cdat[1][0]+cdat[1][1]);
    }else
    {
        nmin2=0.0;
    }

    if(  (MInf>log(numberOfSamples)) && (cdat[0][1]/(cdat[0][0]+cdat[0][1])<=nmin2)  )    {
        return true;
    }else
    {
        return false;
    }
}

/**
 * Initializes the subtree structure and potential sets parents for each gene.
 * (refer to supplementary data section 2)
 *
 * @param DAG       the bayesan network
 * 					Will be modified adding probabilities informations (see initNumberOfMutationProbability)
 * @param tree      the phylogenetic tree initialized with internal nodes and edges.
 * 					Will be modified removing the potential parents the fail the likelihood test.
 */
void localPruning( vector<FAM>& DAG, vector<NODE>& tree )
{

    initializeDAG(DAG);


    for(unsigned int i=0; i<DAG.size()-1; i++)          // for each gene
    {
        if(DAG[i].hasParent==true)                      // for each potential parent
        {
            for(unsigned int j=i+1; j<DAG.size(); j++)  // check following gene
            {
                if( DAG[i].potentialParents[j]==true)   // if i can potentially be a parent of j
                {
                    if(likelihoodTestLocal(i,j,tree)==true)       // and they pass the likelihood test
                    {
                        DAG[j].potentialParents[i]=true;          // j can also be a parent of i
                    }
                    else
                    {
                        DAG[j].potentialParents[i]=false;         // else it can't be
                    }
                }
            }
        }
    }


    initNumberOfMutationProbability(DAG,tree);

}

/**
 * creates a copy in DAG of BN
 * @param DAG     destination
 * @param BN      source
 */
void initDAG(vector<FAM>& DAG, vector<FAM>& BN)
{
    for(unsigned int i=0;i<BN.size();i++)
    {
        DAG[i]=BN[i];
    }
}

/**
 * Computes likelihood score and compares it to logarithm of the total number of possible pairings
 * (refer to chapter 2.2 of the supplementary informations)
 * (global pruning)
 * (A and B are genes)
 * @param numberOfCoOccurrences     number of times A and B are both mutated in a sample
 * @param Atot                      A's number of mutations
 * @param Btot                      B's number of mutations
 * @param numberOfSamples           total number of samples
 * @return the test result         (true if score>log(totalPairings))
 */
bool likelihoodTestGlobal(double& numberOfCoOccurrences, double& Atot, double& Btot, double& numberOfSamples)
{
    double MinfMax=0.0;
    
    /* Note: A1B1 is the number of occurrences where both A and B  mutated, similar for the other variables */

    double A1B0= (double) 1.0*(Atot-1.0*numberOfCoOccurrences);
    double A0B1= (double) 1.0*(Btot-1.0*numberOfCoOccurrences);
    double A1B1= numberOfCoOccurrences;
    double TotS= (double) 2.0*numberOfSamples -2.0;                               // number of possible pairings
    double A0B0= (double) 1.0*(TotS-A1B0-A0B1-A1B1);
    double norm= A0B0;

    for(int i=0; i<numberOfCoOccurrences-1; i++)
    {
        A0B0=norm-i;
        A1B1=numberOfCoOccurrences+i;
        double MinfA=0.0;
        if(A0B0>0)
        {
            MinfA+= (double)A0B0*log((1.0*TotS*A0B0)/(1.0*(A0B0+A0B1)*(A0B0+A1B0)));
        }
        
        if(A1B0>0)
        {
            MinfA+= (double)A1B0*log((1.0*TotS*A1B0)/(1.0*(A0B0+A1B0)*(A1B1+A1B0)));
            
        }
        if(A0B1>0)
        {
            MinfA+= (double)A0B1*log((1.0*TotS*A0B1)/(1.0*(A0B0+A0B1)*(A1B1+A0B1)));
        }
        if(A1B1>0)
        {
            MinfA+= (double)A1B1*log((1.0*TotS*A1B1)/(1.0*(A1B1+A1B0)*(A1B1+A0B1)));
        }
        if(MinfA>MinfMax){
        	MinfMax=MinfA;
        }
    }
    
    if(MinfMax>log(TotS) )           // if the score is over the logarithm of the total number of possible pairings the test is passed
    {
        return true;
    }else
    {
    	return false;
    }
}

/**
 * Implements the global pruning, as defined in chapter 2 of the supplementary notes,
 * to reduce the search space.
 * The BayesNet is set up to only consider the possible relations accepted by this heuristic.
 *
 * @param BayesNet  used for output, content is replaced
 * @param Tree      a tree with initialized vertexes
 * @param op        if 'y' prints a summary of the operations executed on the standard output and on the output file (Pruning)
 * @param job       the name of the job (used for output)
 */
void globalPruning( vector<FAM>& BayesNet, vector<NODE>& Tree, char op, string& job)
{
    FAM famtemp;					// Bayesan Network for a family

    vector<int> TempCO;             // temp for co-occurring mutations
    vector<vector<int> > CoOccur;   // co-occurring mutations (counts the number of times the genes i,k have the same mutation)
    
    /* initialization of FAM data structure (a Markov blanket) */

    for(unsigned int i=0;i<Tree[0].genotype.size(); i++)
    {
        TempCO.push_back(0);
        famtemp.parents.push_back(false);
    }

    famtemp.potentialParents=famtemp.parents;
    famtemp.hasParent=false;
    famtemp.numberOfParents=0;
    famtemp.numberOfMutations=0.0;
    famtemp.score=0.0;

    /* Bayesan network initialization */

    for(unsigned int i=0;i<Tree[0].genotype.size();i++)
    {
        BayesNet.push_back(famtemp);
        CoOccur.push_back(TempCO);
    }
    
    /* computation of Co-Occurrences matrix and of the freq array */

    vector<double> freq(BayesNet.size());          // number of times a certain gene mutated in samples
    vector<int> Order;

    for(unsigned int i=0;i<BayesNet.size();i++)    // for each gene
    {
        Order.push_back(i);                        // keep the input order by now
        freq[i]=0.0;
        for(unsigned int j=0; j<Tree.size(); j++)  // for each sample
        {
            if(Tree[j].genotype[i]==true){
            	freq[i]+=1.0;                      // if there is a mutation count it
            }
            for(unsigned int k=i+1;k<BayesNet.size();k++)                      // for each following gene
            {
                if((Tree[j].genotype[i]==true) && (Tree[j].genotype[k]==true)) // check if there is a co-occurence of the mutation
                {
                    CoOccur[i][k]++;										   // if there is, count it
                    CoOccur[k][i]++;
                }
            }
            
        }
    }

    /* Sorting by decreasing freq (via BubbleSort)*/

    int Count=0;
    for(unsigned int i=0;i<BayesNet.size()-1;i++)
    {
        for(unsigned int j=i+1;j<BayesNet.size();j++)
        {
            if(freq[i]<freq[j])
            {
                int temp=Order[i];      // swap
                Order[i]=Order[j];
                Order[j]=temp;

                double ftemp=freq[i];
                freq[i]=freq[j];
                freq[j]=ftemp;
            } 
        }
    }

    for(unsigned int i=0;i<BayesNet.size();i++) // Saving the new order
    {
        BayesNet[Order[i]].order=i;
    }
    
    /* FAM nodes initialization */
    
    for(unsigned int i=0;i<BayesNet.size()-1;i++)
    {
        
        for(unsigned int j=i+1;j<BayesNet.size();j++)
        {
            double numberOfCoOccurrences = (double) CoOccur[i][j];
            double numberOfSamples = (double) Tree.size();
            int a=BayesNet[i].order;
            int b=BayesNet[j].order;
            
            /* checks if freq[a] and freq[b] together have a sufficient likelihood score */

            if( likelihoodTestGlobal(numberOfCoOccurrences,freq[a],freq[b],numberOfSamples)==true )
            {
                BayesNet[i].potentialParents[j]=true;       // if the test is passed
                BayesNet[j].potentialParents[i]=true;       // this two genes (nodes) can possibly be related
                BayesNet[i].hasParent=true;
                BayesNet[j].hasParent=true;
                Count++;
            }
            
        }
    }

    /* if op == 'y' this procedure prints a summary of this operation on both the STDOUT and the output file */

    if(op=='y')
    {
        string prn= "output/Pruning"+job;
        fstream prnFile;
        prnFile.open(prn.c_str(),ios::out);
        
        cout<<"Total Num of Edges "<<BayesNet.size()*(BayesNet.size()-1)/2.0<<endl<<"Num of unpruned edges "<<Count<<endl;
        prnFile<<"Total Num of Edges "<<BayesNet.size()*(BayesNet.size()-1)/2.0<<endl<<"Num of unpruned edges"<<Count<<endl;
        
        Count=0;
        for(unsigned int i=0;i<BayesNet.size();i++)  // count number of genes with no parent after global pruning
        {   
            if(BayesNet[i].hasParent==false){
                Count++;
            }
        }

        cout<<"Total number of genes "<<BayesNet.size()<<endl<<" Number of genes with no parent after global pruning "<<Count<<endl;
        prnFile<<"Total number of genes "<<BayesNet.size()<<endl<<"Number of genes with no parent after global pruning "<<Count<<endl;
        prnFile.close();
    }
    
    initNumberOfMutationProbability( BayesNet, Tree );

}

