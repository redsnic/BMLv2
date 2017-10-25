/*
 *  SearchDAGs.cpp
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

/*	Performs local search operations on "DAG" given "Tree" for structure learning.				*/

#include "SearchDAGs.h"

using namespace std;

/**
 * dfs used in k2Search
 *
 * @param conf        color gray
 * @param seen        color black
 * @param label
 * @param cdata
 * @param insta
 */
void dfs(vector<bool>& conf, vector<bool>& seen, int& label,
		vector<vector<double> >& cdata, bool& insta) {
	for (unsigned int nb = 0; nb < conf.size(); nb++) {
		if (conf[nb] == false) {
			conf[nb] = true;
			label += int(pow(2.0, nb));
			double nmin1, nmin2;
			if (cdata[0][label] > 0.0) {
				nmin1 = cdata[0][label] / (cdata[0][label] + cdata[1][label]);
			} else {
				nmin1 = 0.0;
			}
			if (cdata[0][label - int(pow(2.0, nb))] > 0.0) {
				nmin2 = cdata[0][label - int(pow(2.0, nb))]
						/ (cdata[0][label - int(pow(2.0, nb))]
								+ cdata[1][label - int(pow(2.0, nb))]);
			} else {
				nmin2 = 0.0;
			}

			if (nmin1 <= nmin2 && insta == true) {
				if (seen[label] == false) {
					dfs(conf, seen, label, cdata, insta);
				}

			} else {
				insta = false;
			}
			seen[label] = true;
			conf[nb] = false;
			label -= int(pow(2.0, nb));
		}

	}
}



/**
 * Greedy K2 heuristic:
 * given a phylogeny it updated the parent of the node relative to a gene.
 * cf. Cooper & Herskovits (1992).
 *
 * @param FamG  a node of the bayesan network
 * @param DAG   bayesan network
 * @param Tree  Phylogeny
 * @param gene  the gene in analysis
 * @param op    (optimization?)
 */
void k2Search( FAM& FamilyOfGene, vector<FAM>& DAG, vector<NODE>& Tree, int gene, char op)
{
	/* Initialization */

    for(unsigned int j=0;j<DAG.size();j++)
    {
        FamilyOfGene.parents[j]=false;
    }
    FamilyOfGene.numberOfParents=0;
    double theta = FamilyOfGene.numberOfMutations/Tree.size();
    FamilyOfGene.k2Parameters.erase(FamilyOfGene.k2Parameters.begin(),FamilyOfGene.k2Parameters.end());
    
    FamilyOfGene.k2Parameters.push_back(1.0 -theta);
    FamilyOfGene.k2Parameters.push_back(theta);
    FamilyOfGene.score=0.0;


	bool repeat=true;
	while(repeat==true)
	{
		double score=FamilyOfGene.score;
		int add=0;
        double numberOfVertexes=(double) Tree.size();
        
		int numberOfParents=FamilyOfGene.numberOfParents;
        vector<double> famparam(int(pow(2.0,numberOfParents+2))), cst(2);
		repeat=false;
		for(unsigned int i=0;i<DAG.size();i++)   // for each node (gene)
		{
            vector< vector<double> > cdata(2);   // a 2*numberOfParents matrix

			if((DAG[gene].potentialParents[i]==true && DAG[i].potentialParents[gene]==true)|| op=='n')
            {
                if( DAG[gene].order>DAG[i].order && FamilyOfGene.parents[i]==false )  // follow the order
                {	
                    for(int m=0;m<int(pow(2.0,numberOfParents+1));m++)   // for all possible instantiations of the node's parents
                    {
                        cdata[0].push_back(0);                           // zero initialization
                        cdata[1].push_back(0);
                    }
                    cst[0]=0;
                    cst[1]=0;

                    numberOfVertexes=0.0;

                    FamilyOfGene.parents[i]=true;                        // simulate i as parent

                    for(unsigned int j=0;j<Tree.size();j++)              // for each sample
                    {
                        int geneState = Tree[j].genotype[gene];
                        int parentsState=0;
                        int parentCount=0;
                        for(unsigned int k=0;k<Tree[j].genotype.size();k++)   // for each gene
                        {
                            if(FamilyOfGene.parents[k]==true)                 // if it is a parent of the current node
                            {
                                parentsState+=int(pow(2.0,parentCount)) * Tree[j].genotype[k];  // if there is a mutation prepare the index to point
                                parentCount++;                                                  // to a given state of parents (that is unique by the use of pow)
                            }
                        }
                        cdata[geneState][parentsState]++;
                        cst[geneState]++;                                                       // total number of mutations, and non mutations
                        numberOfVertexes++;
                    }

                    FamilyOfGene.parents[i]=false;

                    double Minf=0.0;
                    vector<bool> seen;  // has been seen by the dfs visit
                    for(int l=0;l<int(pow(2.0,numberOfParents+1));l++)     // computation of Fam(C, Pi_{C})
                    {
                        if(cdata[0][l]>0.0){Minf+=cdata[0][l]*log((numberOfVertexes*cdata[0][l])/((cdata[0][l]+cdata[1][l])*(cst[0])));}
                        if(cdata[1][l]>0.0){Minf+=cdata[1][l]*log((numberOfVertexes*cdata[1][l])/((cdata[0][l]+cdata[1][l])*(cst[1])));}
                        seen.push_back(false);
                    }
                    vector<bool> conf;
                    for(int l=0;l<numberOfParents+1;l++)
                    {
                        conf.push_back(false);
                    }


                    int label=0;
                    bool insta=true;

                    if(op=='y'){dfs(conf,seen,label, cdata, insta);}  // if not 'y' skip DFS

                    if( insta==true && Minf-pow(2.0,numberOfParents)*log(numberOfVertexes)>score)  // maximize score
                    {   
                        score=Minf-pow(2.0,numberOfParents)*log(numberOfVertexes);
                        add=i;
                        for(int pin=0;pin<int(pow(2.0,numberOfParents+1));pin++)
                        {
                            if(cdata[0][pin]+cdata[1][pin] >0.0)
                            {
                                famparam[pin]=(cdata[0][pin]+ cst[0]/numberOfVertexes)/(cdata[0][pin]+cdata[1][pin]+1.0);
                                famparam[pin+int(pow(2.0,numberOfParents+1))]=(cdata[1][pin] + cst[1]/numberOfVertexes)/((cdata[0][pin]+cdata[1][pin])+1.0);
                            }
                            else
                            {
                                famparam[pin]=1.0;
                                famparam[pin+int(pow(2.0,numberOfParents+1))]=0.0;
                            }
                        }
                    }
                }
            }
            
		}
		if((score>FamilyOfGene.score))
		{
			FamilyOfGene.parents[add]=true;   // node that maximizes the Fam function (see publication)
            FamilyOfGene.score=score;
			FamilyOfGene.numberOfParents++;
            FamilyOfGene.k2Parameters=famparam;
			repeat=true;
		}        
        
	}
    
}


/**
 * Perform local search operations on the bayesan network given the phylogenetic tree for structure learning.
 *
 * @param DAG    the bayesan network (with localPruning)
 * 				 nodes' parents, score and parameters from the K2 algorithm will be added.
 * @param tree   Pylogenetic tree    (with edges)  (unmodified)
 * @param op     mode for k2Search
 */
void searchDAGs( vector<FAM>& DAG, vector<NODE>& Tree, char op )
{
	bool opt=false;

	/* initialization of Ord (order of OBS) */

	vector<int> Ord(DAG.size());
	for(unsigned int i=0;i<DAG.size();i++)
	{
        Ord[DAG[i].order]=i;                               // sorted by number of mutations
	}

	/* for each gene, do a k2Search (to find the parents) */

    for(unsigned int i=0;i<DAG.size();i++)
	{
        if(DAG[i].hasParent==true || op=='n'){ k2Search(DAG[i],DAG,Tree,i,op); }
	}

	while(opt==false)                                      // search for optimum ordering (by perturbating the bayesan network)
	{
		opt=true;
		for(unsigned int i=0;i<DAG.size()-1;i++)           // for each gene
		{
			if((DAG[Ord[i]].potentialParents[Ord[i+1]]==true) && (DAG[Ord[i+1]].potentialParents[Ord[i]]==true))  // if it and its successor can be potentially parent of each other
            {

                int tempord=Ord[i];
                DAG[Ord[i]].order+=1;
                DAG[Ord[i+1]].order-=1;
                Ord[i]=Ord[i+1];
                Ord[i+1]=tempord;                          // swap the orders
                FAM fg1, fg2;
                fg1=DAG[Ord[i]];
                fg2=DAG[Ord[i+1]];
                
                if(DAG[Ord[i]].parents[Ord[i+1]]==true)    // if (in the given order) i+1 is parent of i
                {
                    k2Search(fg1,DAG,Tree,Ord[i],op);
                }
                k2Search(fg2,DAG,Tree,Ord[i+1],op);
                if((fg1.score+fg2.score)>(DAG[Ord[i]].score+DAG[Ord[i+1]].score))   // check if score increased
                {
                    DAG[Ord[i]]=fg1;
                    DAG[Ord[i+1]]=fg2;                                              // if yes keep the new ordering
                    opt=false;                                                      // optimum still has to be found
                }else
                {
                    int tempord=Ord[i];                                             // else swap back
                    DAG[Ord[i]].order+=1;
                    DAG[Ord[i+1]].order-=1;
                    Ord[i]=Ord[i+1];
                    Ord[i+1]=tempord;
                }
            }
		}
	}
    	
}
