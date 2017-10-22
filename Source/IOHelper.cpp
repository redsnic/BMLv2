/*
 *  IOHelper.cpp
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

/*	I/O Stuff & exception handler										*/

#include "IOHelper.h"

using namespace std; 

/**
 * Initializes the given parameters vectors from SimDAG
 *
 * @param netrel1       Edges present in Tree+OBSDAG
 * @param netrel2       Edges present in OBSDAG
 * @param netrel3       Edges present in SimDAG
 * @param EdList        Edges of the SimDAG, with name of associated gene
 * @param SimDAG
 */
void initNetRel(vector<vector<int> >& netrel1,vector<vector<int> >& netrel2,vector<vector<int> >& netrel3, vector< vector<string> >& edgeList, vector<SimVar>& SimDAG)
{
    for(unsigned int i=0;i<SimDAG.size();i++)
    {
        vector<int> temprel;
        for(unsigned int j=0;j<SimDAG.size();j++)
        {
            temprel.push_back(0);
        }

        netrel1.push_back(temprel);
        netrel2.push_back(temprel);
        netrel3.push_back(temprel);
        
        string name = SimDAG[i].name;
        
        for(unsigned int j=0;j<SimDAG[i].parents.size();j++)
        {
            int parentIndex=SimDAG[i].parents[j];
            netrel3[i][parentIndex]=1;
            netrel3[parentIndex][i]=1;
            vector<string> parentName;
            parentName.push_back(SimDAG[parentIndex].name);
            parentName.push_back(name);
            edgeList.push_back(parentName);
        }
    }
}

/**
 * Output procedure: prints the TruePositive / FalsePositives confidence intervals
 *
 * @param NetRel1
 * @param NetRel2
 * @param NetRel3
 * @param SimDAG
 * @param NRep
 * @param job
 */
void tpfpAnalysis( vector<vector<int> >& NetRel1,vector<vector<int> >& NetRel2, vector<vector<int> >& NetRel3, vector<SimVar>& SimDAG, int& NRep, string& job)
{

	/* Opening output files */

    fstream relFile;
    string tp="output/TruePositives"+job;
    string fp= "output/FalsePositives"+job;
    relFile.open(tp.c_str(),ios::out);
    relFile<<'\t'<<"Tree+DAG "<<'\t'<<"DAG"<<endl;
    vector<int> Hist1, Hist2;        // histograms
    
    fstream TruePositivesConfidence;
    string tpc="output/ConfidenceTruePositives"+job;
    TruePositivesConfidence.open(tpc.c_str(),ios::out);
    
    fstream FalsePositivesConfidence;
    string fpc="output/ConfidenceFalsePositives"+job;
    FalsePositivesConfidence.open(fpc.c_str(),ios::out);

    /* histogram creation */
    
    for(int i=0;i<10;i++)                              // initialization // TODO use a constant instead of 10
    {
        Hist1.push_back(0); 
        Hist2.push_back(0); 
    }

    /* histogram */

    for(unsigned int i=0;i<NetRel1.size()-1;i++)
    {
        for(unsigned int j=i+1;j<NetRel1.size();j++)   // population
        {
            if(NetRel3[i][j]==true)                    // for each edge
            {
                relFile<<SimDAG[i].name<<"-"<<SimDAG[j].name<<'\t';
                relFile<<NetRel1[i][j]<<'\t'<<NetRel2[i][j];
                relFile<<endl;
                int bin1=int(10*NetRel1[i][j]/NRep);
                int bin2=int(10*NetRel2[i][j]/NRep);
                if(bin1==10){bin1=9;}
                if(bin2==10){bin2=9;}
                Hist1[bin1]++;
                Hist2[bin2]++;
            }
        }
    }

    /* true positives */

    relFile<<endl<<endl;
    relFile.close();
    TruePositivesConfidence<<"Confidence"<<'\t'<<"Tree+DAG"<<'\t'<<"DAG"<<endl;
    int nted=0;
    int noed=0;
    for(int i=9;i>0;i--)  // TODO constant
    {
        nted+=Hist1[i];
        noed+=Hist2[i];
        TruePositivesConfidence<<10*i<<'\t'<<nted<<'\t'<<noed<<endl;
        Hist1[i]=0; 
        Hist2[i]=0; 
    }
    Hist1[0]=0;
    Hist2[0]=0;
    TruePositivesConfidence.close();
    
    fstream relFile2;
    relFile2.open(fp.c_str(),ios::out);
    relFile2<<'\t'<<"Tree+DAG "<<'\t'<<"DAG"<<endl;
    
    /* histogram */

    for(unsigned int i=0;i<NetRel1.size()-1;i++)
    {
        
        for(unsigned int j=i+1;j<NetRel1.size();j++)
        {
            if((NetRel1[i][j]+NetRel2[i][j]>0)&&(NetRel3[i][j]==0))   // histogram preparation
            {
                relFile2<<SimDAG[i].name<<"-"<<SimDAG[j].name<<'\t';
                relFile2<<NetRel1[i][j]<<'\t'<<NetRel2[i][j];
                relFile2<<endl;
                int bin1=int(10*NetRel1[i][j]/NRep);
                int bin2=int(10*NetRel2[i][j]/NRep);
                if(bin1==10){bin1=9;}
                if(bin2==10){bin2=9;}
                Hist1[bin1]++;
                Hist2[bin2]++;
            }
        }
    }
    relFile2.close();
    
    /* false positives */

    nted=0;
    noed=0;
    FalsePositivesConfidence<<"Confidence"<<'\t'<<"Tree+DAG"<<'\t'<<"DAG"<<endl;
    for(int i=9;i>0;i--)
    {
        nted+=Hist1[i];
        noed+=Hist2[i];
        FalsePositivesConfidence<<10*i<<'\t'<<nted<<'\t'<<noed<<endl;
        Hist1[i]=0; 
        Hist2[i]=0; 
    }
    Hist1[0]=0;
    Hist2[0]=0;
    FalsePositivesConfidence.close();
    
}

/**
 * Computes the probability ratio of a parent configuration
 *
 * @param SimDAG               output mutations DAG
 * @param config               config array ( the probability is computed for the set of true elements )
 * @param probabilityRatio     output (this is the only modified parameter)
 */
void configParam(vector<SimVar>& SimDAG, vector<bool>& config, double& probabilityRatio)
{
    for(unsigned int l=0;l<SimDAG.size();l++)
    {
        int y=0;
        int y2=SimDAG[l].parents.size();
        for(int pa=0;pa<y2;pa++)
        {
            y+=config[SimDAG[l].parents[pa]]*(int(pow(2.0,pa)));
        }
        double z=SimDAG[l].k2Parameters[y]/(SimDAG[l].k2Parameters[y]+SimDAG[l].k2Parameters[y+(int(pow(2.0,y2)))] );
        if(config[l]==true)
        {
        	probabilityRatio=probabilityRatio*(1.0-z);
        }
        else{
        	probabilityRatio=probabilityRatio*z;
        }
    }
}


/**
 * output procedure: maximum depth of DAG is set to 3
 * The procedure gets the triplets of genes'mutations with a probability (normalized to the maximum probability of a triplet)
 * greater than the threshold set by the user. Then it reconstructs the most probable parents and does the same for singolet and pairs of genes.
 *
 * @param SimDAG         DAG for output
 * @param tree           Phylogeny
 * @param GeneLabels     Genes' name
 * @param char_label     Genes' groups
 * @param pthres         threshold given as input
 * @param job            job name for output
 */
void writeLandscape(vector<SimVar>& SimDAG, vector<NODE>& tree, double& pthres, string& job)
{
    const double epsilon=0.0001;

    int NVars= SimDAG.size();
	vector<double> NormStat;
	vector<bool> config(SimDAG.size());
    
	/* opening output file */

    fstream pathFile;
    string bpath="output/Path"+job+".dot";
	pathFile.open(bpath.c_str(),ios::out);


	for(unsigned int i=0;i<SimDAG.size();i++)
	{
		config[i]=false;
	}

	for(unsigned int i=0;i<SimDAG.size();i++)   // for each gene
	{
		double tempsc=0.0;
		tempsc=SimDAG[i].k2Parameters[0];       // take first parameter
		NormStat.push_back(tempsc);             // and store it for normalization purpose
	}
    
	/* preparing .dot file */
    
    pathFile << "digraph G{" << endl << "rankdir= LR;" << endl; // header
    
    pathFile<<0<<"[label=\"Normal\" ,shape=ellipse, style=filled, color=\"0.000,1.0,0.85\",fontsize=\"36\", font=\"Helvetica\", fontcolor=white];"<<endl;
    
    double max=0.0;
    vector<double>  prob;
    vector<string> plabel;
    vector<bool> present;
	for(unsigned int i=0;i<SimDAG.size();i++)               // for each gene
	{
		config[i]=true;
		double probabilityRatio=1.0;
        configParam( SimDAG,  config,  probabilityRatio);

        if(probabilityRatio>max)
		{
			max=probabilityRatio;
		}
        prob.push_back(probabilityRatio);                   // compute its probability ratio

        present.push_back(false);
		config[i]=false;
	}
    
    vector<double> prob2;
    vector< vector<int> > plabel2;
    vector<bool> present2;
    double max2=0.0;
	for(unsigned int i=0;i<SimDAG.size()-1;i++)             // for each pair of genes (i<j)
	{
		config[i]=true;
        for(unsigned int j=i+1;j<SimDAG.size();j++)
        {	
            config[j]=true;
            
            double pratio=1.0;		
            configParam( SimDAG,  config,  pratio);        // compute their probability ratio
            
            
            prob2.push_back(pratio);
            vector<int> lab;
            lab.push_back(i);
            lab.push_back(j);
            plabel2.push_back(lab);
            present2.push_back(false);
            
            if(pratio>max2){max2=pratio;}
            
            config[j]=false;
        }
		config[i]=false;
	}

    vector<double> prob3;                                   // for each triplet of genes (i<j<k)
    vector< vector<int> > plabel3;
    vector<bool> present3;
    double max3=0.0;
	for(unsigned int i=0;i<SimDAG.size()-2;i++)
	{
		config[i]=true;
        for(unsigned int j=i+1;j<SimDAG.size()-1;j++)
        {	
            
            config[j]=true;
            
            
            for(unsigned int k=j+1;k<SimDAG.size();k++)
            {
                config[k]=true;
                bool pre=false;
                for(int pr=1;pr<tree[0].numberOfDescendantSamples+1;pr++)
                {
                    if(tree[pr].genotype[SimDAG[i].geneIndex]==true && tree[pr].genotype[SimDAG[j].geneIndex]==true && tree[pr].genotype[SimDAG[k].geneIndex]==true)
                    {
                        pre=true;
                    }
                }
                present3.push_back(pre);                    // if the mutation is present in all the three nodes push true in present3
                double pratio=1.0;		
                configParam( SimDAG,  config,  pratio);     // compute their probability ratio
                
                
                prob3.push_back(pratio);
                vector<int> label;
                label.push_back(i);
                label.push_back(j);
                label.push_back(k);
                plabel3.push_back(label);                    // label for the new triplet node
                
                if(pratio>max3){max3=pratio;}
                config[k]=false;
                
            }
            config[j]=false;
        }
		config[i]=false;
	}
    for(unsigned int i=0;i<prob3.size();i++)
    {
        
        if( present3[i]==true && prob3[i]>=max3*(pthres-epsilon)) // takes the pthres% more probable mutations
        {
            pathFile<<i+1+NVars+NVars*(NVars-1)/2<<"[label=\""<<SimDAG[plabel3[i][0]].name<<","<<SimDAG[plabel3[i][1]].name<<","<<SimDAG[plabel3[i][2]].name<< " : " << prob3[i]/max3 <<"\" ,shape=ellipse, style=filled, color=\"0.13000,"<<prob3
            [i]/max3<<",0.99\",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];"<<endl;
            int pa = 0;
            double papr=0;
            for(unsigned int j=0; j<prob2.size(); j++)
            {
                if( ( (plabel2[j][0]==plabel3[i][0]) && ( (plabel2[j][1]==plabel3[i][1]) || (plabel2[j][1]==plabel3[i][2]) ) )  // two labels must be the same of a parent
                		|| ( (plabel2[j][0]==plabel3[i][1]) &&  (plabel2[j][1]==plabel3[i][2])  ) )
                {
                    if(prob2[j]>papr)
                    {
                        pa=j;
                        papr=prob2[j];
                    }
                }
            }
            present2[pa]=true;
            pathFile<<NVars+pa+1<<"->"<<i+1+NVars+NVars*(NVars-1)/2<<endl;  // prints edge
            
        }
    }
    
    for(unsigned int i=0;i<prob2.size();i++)
    {
        int g1 = plabel2[i][0];
        int g2 = plabel2[i][1];
        if(prob2[i]>=max2*(pthres-epsilon))
        {
            for(int pr=1;pr<tree[0].numberOfDescendantSamples+1;pr++)
            {
                bool pr2=false;
                
                if(tree[pr].genotype[SimDAG[g1].geneIndex]==true && tree[pr].genotype[SimDAG[g2].geneIndex]==true )   // if they both have a mutation
                {
                    pr2=true;
                }
                int mutationCount=0;
                for(unsigned int j=0;j<tree[pr].genotype.size();j++)
                {
                    if( tree[pr].genotype[j]==true){mutationCount++;}      // and the number of mutations is 2
                }
                if(pr2 == true && mutationCount == 2){present2[i]=true;}   // add edges
            }
        }
        
        if(present2[i]==true )
        { 
            pathFile<<NVars+i+1<<"[label=\""<<SimDAG[g1].name<<","<<SimDAG[g2].name<< " : " << prob2[i]/max2 <<"\" ,shape=ellipse, style=filled, color=\"0.11000,"<<prob2[i]/max2<<",0.90\",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];"<<endl;
            if(prob[g1]>prob[g2])
            {
                pathFile<<g1+1<<"->"<<NVars+i+1<<endl;
                present[g1]=true;
            }else
            {
                pathFile<<g2+1<<"->"<<NVars+i+1<<endl;
                present[g2]=true;
                
            }
        }
        
    }
    
    for(unsigned int i=0;i<prob.size();i++)   // single mutation
    {
        if(prob[i]>=max*(pthres-epsilon))
        {
            for(int pr=1;pr<tree[0].numberOfDescendantSamples+1;pr++)
            {
                bool pr1=false;
                
                if(tree[pr].genotype[SimDAG[i].geneIndex]==true )
                {
                    pr1=true;
                }
                int mutationCount=0;
                for(unsigned int j=0;j<tree[pr].genotype.size();j++)
                {
                    if( tree[pr].genotype[j]==true){mutationCount++;}
                }
                if(pr1==true && mutationCount ==1){present[i]=true;}
            }
        }
        
        if(present[i]==true)
        {
            pathFile<<i+1<<"[label=\""<<SimDAG[i].name<< " : " << prob[i]/max <<"\" ,shape=ellipse, style=filled, color=\"0.075000,"<<prob[i]/max <<",0.90\",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];"<<endl;
            pathFile<<0<<"->"<<i+1<<endl;
        }
        
    }
    
    /* close graph and file */

    pathFile<<"}";
	pathFile.close();
}

/**
 * Computes the probability of each edge in the bayesan network
 *
 * @param SimDAG        output DAG
 * @param landrel       set of probabilities for each edge of the bayesan network
 * @param DAG           bayesan network
 */
void landBoot(vector<SimVar>& SimDAG, vector< vector<double> >& landrel, vector<FAM>& DAG)
{
    int nodeNumber=SimDAG.size();
    vector<string>VarName(nodeNumber);
    vector<vector<int> > VarPi(nodeNumber);
    vector<vector<double> >VarParam(nodeNumber);
	vector<bool>config(DAG.size());

	/* Variables initialization */

    for(int i=0;i<nodeNumber;i++)                     // for each node
    {
    	string vname;
        int ord,numberOfParents;
        vname = SimDAG[i].name;
        ord=DAG[i].order;
        numberOfParents=DAG[i].numberOfParents;
        VarName[ord]=vname;

        for(int j=0;j<nodeNumber;j++)
        {
            if(DAG[i].parents[j]==true){              // for each node
                VarPi[ord].push_back(DAG[j].order);   // store its parents
            }
        }

        int numberOfParameters=int(pow(2.0,numberOfParents+1));
        for(int j=0;j<numberOfParameters;j++)
        {
            double vparam=DAG[i].k2Parameters[j];       // get k2Parameters
            VarParam[ord].push_back(vparam);
        }
    }

    double pnor=1.0;
	for(unsigned int i=0;i<DAG.size();i++)
	{
		config[i]=false;
        pnor=pnor*DAG[i].k2Parameters[0];               // compute probability of a network without any relation
	}

    vector<double>  prob;
    vector<string> plabel;
	for(int ng=0;ng<nodeNumber;ng++)                    // for each node
    {
        for(unsigned int i=0;i<DAG.size();i++)          // for each gene
        {
            if(VarName[i]==SimDAG[ng].name)             // for the gene in the given order
            {
                config[i]=true;
                double pratio=1.0;
                for(unsigned int l=0;l<DAG.size();l++)  // for each gene
                {
                    int parentState=0;
                    int numberOfParents=VarPi[l].size();

                    for(unsigned int pa=0;pa<VarPi[l].size();pa++)   // for each parent
                    {
                        parentState+=config[VarPi[l][pa]]*(int(pow(2.0,pa)));  // encode parent state
                    }

                    double z=VarParam[l][parentState]/(VarParam[l][parentState]+VarParam[l][parentState+(int(pow(2.0,numberOfParents)))] );

                    if(config[l]==true){
                    	pratio=pratio*(1.0-z);          // 1-p
                    }
                    else{
                    	pratio=pratio*z;                // p
                    }
                }

                prob.push_back(pratio);
                config[i]=false;
            }
        }
    }
    
    for(unsigned int i=0;i<prob.size();i++)             // probability for each edge
    {
        landrel[i].push_back(prob[i]/(pnor));
    }
    
}

/**
 * Computes the Bootstrap probabilities for each pair of genes sharing an edge
 * for the 2 single mutant states and the double mutant genotype using Tree+OBSDAG
 *
 * @param edgeList   labelled edges
 * @param SimDAG     the output DAG
 * @param edrel      output, edge probabilities
 * @param DAG		 the bayesan network
 */
void edgeBoot(vector< vector<string> >& edgeList, vector<SimVar>& SimDAG, vector< vector<double> >& edrel, vector<FAM>& DAG)
{

	/* variables initialization */

    int NVars=DAG.size();
    vector<string>VarName(NVars);
    vector<vector<int> > VarPi(NVars);
    vector<vector<double> >VarParam(NVars);
	vector<bool>config(DAG.size());


    for(int i=0;i<NVars;i++)
    {
        string vname;
        int ord,npi;
        vname=SimDAG[i].name;
        ord=DAG[i].order;
        npi=DAG[i].numberOfParents;
        VarName[ord]=vname;
        for(int j=0;j<NVars;j++)
        {
            if(DAG[i].parents[j]==true){
                VarPi[ord].push_back(DAG[j].order);
            }
        }
        int nparam=int(pow(2.0,npi+1));
        for(int j=0;j<nparam;j++)
        {
            double vparam=DAG[i].k2Parameters[j];
            VarParam[ord].push_back(vparam);
        }
    }

    /* probability computation */

    double pnor=1.0;

    for(unsigned int i=0;i<DAG.size();i++)
	{
		config[i]=false;
        pnor=pnor*DAG[i].k2Parameters[0];
	}

	for(unsigned int ng=0;ng<edgeList.size();ng++)      // for each edge
    {
        vector<double> prob(3);
        
        vector<int> glab(2);
        
        for(unsigned int i=0;i<DAG.size();i++)          // for each gene
        {
            if(VarName[i]==edgeList[ng][0] || VarName[i]==edgeList[ng][1] )   // that is adjacent to the edge
            {
                
                config[i]=true;                                               // true for the elements adjacent to the edge
                double pratio=1.0;

                for(unsigned int l=0;l<DAG.size();l++)   // for each other gene
                {
                    int y=0;
                    int numberOfParents=VarPi[l].size();
                    for(unsigned int pa=0;pa<VarPi[l].size();pa++)
                    {
                        y+=config[VarPi[l][pa]]*(int(pow(2.0,pa)));
                    }
                    
                    double z=VarParam[l][y]/(VarParam[l][y]+VarParam[l][y+(int(pow(2.0,numberOfParents)))]);
                    
                    if(config[l]==true){
                    	pratio=pratio*(1.0-z);       // 1-p
                    }
                    else{
                    	pratio=pratio*z;             // p
                    }
                }

                if(VarName[i]==edgeList[ng][0])      // if it is the source
                {
                    glab[0]=i;
                    prob[0]=pratio;                  // first mutated probability
                }else                                // if it is the destination
                {
                    glab[1]=i;
                    prob[1]=pratio;                  // second mutated probability
                }
                config[i]=false;
            }
        }
        config[glab[0]]=true;
        config[glab[1]]=true;
        
        double pratio=1.0;
        for(unsigned int l=0;l<DAG.size();l++)        // for each gene
        {
            int y=0;
            int numberOfParents=VarPi[l].size();
            
            for(unsigned int pa=0;pa<VarPi[l].size();pa++)
            {
                y+=config[VarPi[l][pa]]*(int(pow(2.0,pa)));
            }
            double z=VarParam[l][y]/(VarParam[l][y]+VarParam[l][y+(int(pow(2.0,numberOfParents)))]);
            
            if(config[l]==true){
            	pratio=pratio*(1.0-z);
            }
            else{
            	pratio=pratio*z;
            }
        }
        prob[2]=pratio;                              // both mutated probability
        
        config[glab[1]]=false;
        config[glab[0]]=false;
        for(int bp=0;bp<3;bp++)
        {
            edrel[3*ng+bp].push_back((prob[bp]/pnor));
        }
    }
}

/**
 * Output procedure: edge probabilities for the final DAG
 *
 * @param edgeList   list of edges
 * @param edrel      edges' probabilities
 * @param job        output file suffix
 */
void probEdBoot(vector< vector<string> >& edgeList, vector< vector<double> >& edrel, string& job)
{
    string edb= "output/EdgeProbabilities"+job;
    fstream edFile;
    edFile.open(edb.c_str(),ios::out);
    for(unsigned int ng=0;ng<edgeList.size();ng++)
    {
        edFile<<edgeList[ng][0];
        for(unsigned int bp=0;bp< edrel[0].size();bp++)
        {
            edFile<<'\t'<<edrel[3*ng][bp];
        }
        edFile<<endl;
        edFile<<edgeList[ng][1];
        for(unsigned int bp=0;bp< edrel[0].size();bp++)
        {
            edFile<<'\t'<<edrel[3*ng+1][bp];
        }
        edFile<<endl;
        edFile<<edgeList[ng][0]<<","<<edgeList[ng][1];
        for(unsigned int bp=0;bp< edrel[0].size();bp++)
        {
            edFile<<'\t'<<edrel[3*ng+2][bp];
        }
        edFile<<endl;
        edFile<<edgeList[ng][0]<<","<<edgeList[ng][1]<<"(indep)";
        for(unsigned int bp=0;bp< edrel[0].size();bp++)
        {
            edFile<<'\t'<<edrel[3*ng+1][bp]*edrel[3*ng][bp];
        }
        edFile<<endl;
    }
    edFile.close();
}

/**
 * Output procedure: prints the computed edge probabilities for OBS and OBS+Tree
 *
 * @param NetRel1
 * @param NetRel2
 * @param SimDAG
 * @param job        output file suffix
 */
void probBootstrapAnalysis( vector<vector<double> >& NetRel1,vector<vector<double> >& NetRel2, vector<SimVar>& SimDAG, string& job)
{
    string lr1= "output/OBS_Probabilities"+job;
    string lr2= "output/Tree_OBS_Probabilities"+job;
    
    fstream relFile1;
    fstream relFile2;
    
    relFile1.open(lr1.c_str(),ios::out);
    relFile2.open(lr2.c_str(),ios::out);
    
    for(unsigned int i=0;i<NetRel1.size();i++)
    {
        relFile1<<SimDAG[i].name<<'\t';
        relFile2<<SimDAG[i].name<<'\t';

        for(unsigned int j=0;j<NetRel1[0].size();j++)
        {
            relFile1<<NetRel1[i][j]<<'\t';
            relFile2<<NetRel2[i][j]<<'\t';
        }
        relFile1<<endl;
        relFile2<<endl;
    }
    
    relFile1.close();
    relFile2.close();
}


/**
 * Output the local search based MPE estimates for DAG structure for each starting tree in .dot format
 * used only for debug purpose
 *
 * @param DAG
 * @param GeneLabels
 * @param char_label
 * @param op               'y' to output to Tree_DAG, otherwise it outputs to OBS_DAG
 * @param job
 */
void writeDAG(vector<FAM>& DAG, vector<string>& GeneLabels,
        vector<vector<int> >& char_label, char op) {

	/* opening files */

	ofstream DAGFile;
	fstream outFile;

    string st1 = "output/Tree_DAG";
	string dt1 = st1 + ".dot";
    string st2 = "output/OBS_DAG";
	string dt2 = st2 + ".dot";

	if (op == 'y') {
		DAGFile.open(st1.c_str(), ios::out);
		outFile.open(dt1.c_str(), ios::out); //prepare the output dot file

	} else {
		DAGFile.open(st2.c_str(), ios::out);
		outFile.open(dt2.c_str(), ios::out); //prepare the output dot file

	}

	outFile << "digraph G{" << endl;
	DAGFile << DAG.size() << '\n';

	vector<bool> cflag(DAG.size());

	for (unsigned int i = 0; i < cflag.size(); i++) {  // initialize cflag
		cflag[i] = false;
	}

	for (unsigned int i = 0; i < DAG.size(); i++) {                   // informations about edges
		for (unsigned int j = 0; j < DAG[i].parents.size(); j++) {    // for each node draw the edge from the parents to it
			//int wt=0;
			if (DAG[i].parents[j] == true) {
				outFile << j << "->" << i << endl;
				cflag[i] = true;
				cflag[j] = true;
			}

		}

	}
	for (unsigned int i = 0; i < DAG.size(); i++) {         // informations about vertex
		if (cflag[i] == true) {                             // for each node that has a parent set label and style
			outFile << i
					<< "[style=filled, color=\"skyblue2\",font=\"Helvetica\",label=\""
					<< GeneLabels[char_label[i][0]];
			for (unsigned int k = 1; k < char_label[i].size(); k++) {
				outFile << "," << GeneLabels[char_label[i][k]];
			}
			outFile << "\"]" << endl;
		}

		/* print IDs and informations about each node */

		DAGFile << GeneLabels[char_label[i][0]] << '\t' << DAG[i].order << '\t'
				<< DAG[i].numberOfParents << '\t';

		for (unsigned int j = 0; j < DAG[i].parents.size(); j++) {
			if (DAG[i].parents[j] == 1) {
				DAGFile << DAG[j].order << '\t';
			}

		}

		for (unsigned int j = 0; j < DAG[i].k2Parameters.size(); j++) {
			DAGFile << DAG[i].k2Parameters[j] << '\t';

		}
		DAGFile << endl;

	}

	outFile << "}" << endl;
	outFile.close();
	DAGFile.close();

}



/*	Exception handling and error messages								*/

/**
 * Prints informations about a given exception
 *
 * @param ex         an exception (a number among 0,1,2)
 */
void exceptionHandler(int& ex)
{
	switch (ex) 
	{
		case 0:
			cerr<<"Could not open the data file. Perhaps you entered an incorrect path or file name.\nBailing Out !\n";
			break;
		case 1:
			cerr<<"Please check the data matrix file.\nSee README file or look at the examples in the data folder for detailed formatting instructions.\nBailing Out !\n";
			break;
        case 2:
			cerr<<"child-parent relations awry. Failed the sanity check in InitTree.cpp";
			break;
	}	
}
