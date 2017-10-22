/*
 *  SearchTrees.cpp
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

/*	Perform local search operations in tree space for imputing missing data values
 *	in variable Tree using the MPE (Most Probable Explanation) criterion.
 */

#include "SearchTrees.h"

using namespace std; 

/**
 * the problem is solved by adding a mutation for that gene in the node with the lowest number of descendant samples that
 * has both children with a mutation for that gene.
 *
 * @param tree                    Phylogenetic tree (will be updated correcting the error)
 * @param gene                    gene in wrong position
 */
void relabelOutlier(vector<NODE>& tree,  int gene)
{
    int min=tree.size();
    int label=0;
    for(unsigned int i=tree.size()/2+1;i<tree.size();i++)      // for each internal node
    {
        if(tree[i].genotype[gene]==false && (tree[tree[i].children[0]].genotype[gene]*tree[tree[i].children[1]].genotype[gene]==true))  // new mutation on both children
           {
               if(min>tree[i].numberOfDescendantSamples)
               {
                   min=tree[i].numberOfDescendantSamples;
                   label=i;                                    // take the node with the minimum number of descendant samples (respecting the if conditions)
               }
           }
    }
    tree[label].genotype[gene]=true;                           // and set it to have that gene as mutated
}

/**
 * Swap the i node with its uncle preserving the correctness of the phylogenetic tree.
 *
 * @param Tree            phylogenetic tree (nodes will be swapped)
 * @param DAG             bayesan network   (local pruning will be updated accordingly to the new phylogenetic tree)
 * @param i               considered node
 * @param grandparent
 * @param uncle
 * @param change          output: true if there was a change in internal nodes
 */
void swapnodes(vector<NODE>& tree, vector<FAM>& DAG, int& i, int& grandparent, int& uncle, bool& change)
{
    int isNodeARightChild=0;
    int sibling=tree[tree[i].parent].children[1];

    if(tree[tree[i].parent].children[1]==i)                              // if it is a right child, its sibling is the left node
    {
        isNodeARightChild=1;
        sibling=tree[tree[i].parent].children[0];
    }

    int isParentALeftChild=0;
    if(tree[grandparent].children[0]==tree[i].parent){isParentALeftChild=1;}
    
    uncle = tree[grandparent].children[isParentALeftChild];
    
    change=false;
    for(unsigned int j=0;j<tree[i].genotype.size();j++)                                 // for each gene
    {                                                                                   // check mutations in:
        bool uncleAndSiblingState=tree[uncle].genotype[j]*tree[sibling].genotype[j];    // uncle AND sibling
        bool parentState=tree[tree[i].parent].genotype[j];                              // parent
        bool thisAndSiblingState = tree[i].genotype[j]*tree[sibling].genotype[j];       // this sample AND its sibling
        if((thisAndSiblingState!=parentState))                                          // if there are differences with the parent:
        																			    // it is possible that the parent of the two nodes had a mutation and their children had not
        {                                                                               // (which is impossible by hypothesis), or the other way around (both children mutated).
            if((DAG[j].numberOfMutations==tree.size()/2) && parentState==false)         // if both children mutated and the number of mutations corresponds to the number of leaves
            {
                uncleAndSiblingState=false;                                             // parent may be updated eliminating mutation
            }
            else                                                                        // else there is an error in tree creation
            {
                    cout<<"incorrect label for node "<<i<<" sibling "<<sibling<<" dad "<<tree[i].parent<<" grandparent "<<grandparent<<" uncle "<<uncle<<endl;
                    cout<<"Nmut "<<DAG[j].numberOfMutations<<endl;
                    cout<<"Dad's children "<<tree[tree[i].parent].children[0]<<" and "<<tree[tree[i].parent].children[1]<<endl;
                    cout<<"grandparent children "<<tree[grandparent].children[0]<<" and "<<tree[grandparent].children[1]<<endl;
                    cout<<"uncle's father "<<tree[uncle].parent<<endl;
                    cout<<"sibling's father "<<tree[sibling].parent<<endl;
                    cout<<j<<" gene, node ="<<tree[i].genotype[j]<<" sibling = "<<tree[sibling].genotype[j]<<" dad = "<<tree[tree[i].parent].genotype[j]<<" gramps = "<<tree[grandparent].genotype[j]<<" uncle = "<<tree[uncle].genotype[j]<<endl;
                    throw 2;
            }
        }
        double theta=0.0;
        if(parentState!=uncleAndSiblingState)                                             // if there is a change from the parent mutations state in uncle and sibling
        {
            change=true;
            if(parentState==true)                                                         // if there was a mutation in parent
            {
                if(DAG[j].numberOfMutations==tree.size()/2)								  // and if there are yet enough mutations
                {
                    relabelOutlier(tree, j);                                              // correct anomaly (mutations must be permanent)
                }else
                {
                	DAG[j].numberOfMutations-=1.0;                                        // else decrease number of mutations (a mutation will be removed)
                }
            }else                                                                         // if there was not a mutation in parent (so both sibling and uncle are mutated)
            {
                if(DAG[j].numberOfMutations<tree.size()/2)                                // and if the number of mutations is inferior to the number of leaves
                {
                	DAG[j].numberOfMutations+=1.0;										  // increase the number of mutations (a mutation will be added)
                }
                else
                {
                	uncleAndSiblingState=false;                                           // else number of mutation should be kept the same
                }
            }
            theta = DAG[j].numberOfMutations/tree.size();                   // mutation probability for gene j
            if(theta>0.0)
            {
                DAG[j].numberOfMutationsProbability= -DAG[j].numberOfMutations*log(theta) - (tree.size()-DAG[j].numberOfMutations)*log(1.0-theta);   // log probability
            }else
            {
                DAG[j].numberOfMutationsProbability=0.0;                    // avoid sub zero probabilities
            }

            tree[tree[i].parent].genotype[j]=uncleAndSiblingState;          // update parent genotype

            for(unsigned int k=0;k<j; k++)                                  // for each gene before j
            {                                                               // (update local pruning using the new tree)
                if(DAG[k].potentialParents[j]==true)                        // decide which one can still be a potential parent of the other
                {
                    if(likelihoodTestLocal(k, j, tree)==true)
                    {
                        
                            DAG[j].potentialParents[k]=true;
                        }
                        else
                        {
                            DAG[j].potentialParents[k]=false; 
                        }

                    
                }
                if(DAG[j].potentialParents[k]==true)
                {
                    if(likelihoodTestLocal(j, k,tree )==true)
                    {
                        
                            DAG[k].potentialParents[j]=true;
                        }
                        else
                        {
                            DAG[k].potentialParents[j]=false; 
                        }
                    
                    
                }
                
            }
            
        }
    }
    
    tree[grandparent].children[isParentALeftChild]=i;                          // swap the given node with its uncle
    tree[tree[i].parent].children[isNodeARightChild]=uncle;
    tree[uncle].parent=tree[i].parent;
    tree[i].parent=grandparent;
}


/**
 * Searches for valid trees in the search space
 *
 * @param DAG     Boolean Network (with localPruning).
 * 				  Nodes' order and relative parents, scores and parameters will be updated.
 * @param tree    phylogeny tree.
 * 				  It will be modified with the tree locally maximizing the log-likelihood score.
 */
void searchTrees( vector<FAM>& DAG, vector<NODE>& tree)
{

	bool opt =false;
    double oldScore=0.0;
    int degree=0;

    searchDAGs(DAG, tree, 'y');

	for(unsigned int i=0; i<DAG.size();i++)       // for each gene (family, Markov blanket)
	{
		oldScore+=DAG[i].score-DAG[i].numberOfMutationsProbability;       // compute tree score (log-liklehood)
        degree+=DAG[i].numberOfParents;                                   // and degree
	}

	while(opt==false)                             // search for optimal score (by tree perturbation)
	{
		opt=true;
        cout<<" tree size = "<<tree.size()<<" present score "<< oldScore <<'\n';

		for(unsigned int i=1; i<tree.size();i++)  // for each sample
		{

            sanityCheckTree(tree);                // check for errors in tree
            int uncle=0;
            int candidate=i;
            double newScore=0.0;
            vector<int> tempOrder;
			int grandparent = tree[tree[i].parent].parent;
            bool change=false;
			if(grandparent!=0)                                               // if the grandparent is not the the root
			{
				swapnodes(tree, DAG, candidate, grandparent, uncle, change); // swap the candidate with its uncle
				if (change == true)                                          // if nodes were swapped
				{
					for (unsigned int j = 0; j < DAG.size(); j++) {
						tempOrder.push_back(DAG[j].order);                   // save old order
					}
					searchDAGs(DAG, tree, 'y');                              // new research with modified tree

					for (unsigned int j = 0; j < DAG.size(); j++)
					{
						newScore += DAG[j].score - DAG[j].numberOfMutationsProbability; // compute new DAG score
					}
					if (newScore > oldScore)                                            // if the new order is better
					{
						opt = false;			                                        // optimum still has to be found
						oldScore = newScore;			                          		// update score
						degree = 0;
						for (unsigned int j = 0; j < DAG.size(); j++)
						{
							degree += DAG[j].numberOfParents;                           // compute new degree
						}
					} else {
						swapnodes(tree, DAG, uncle, grandparent, candidate, change);    // else swap back
						for (unsigned int j = 0; j < DAG.size(); j++)
						{
							DAG[j].order = tempOrder[j];                                // reset order
						}
					}
				}
			}
		}
	}

    double finalScore=0.0;
        
    /* print operation's results */

	searchDAGs(DAG, tree,'y');
    int nedges=.0;
    for(unsigned int j=0; j<DAG.size();j++)
    {
        finalScore+=DAG[j].score-DAG[j].numberOfMutationsProbability;
        nedges+=DAG[j].numberOfParents;
    }

    cout << nedges << " edges and log-likl = " << finalScore << '\n';


}
