#include "LandscapeTree.h"

LandscapeTree::LandscapeTree(std::vector<SimVar>& DAG, double pthres, std::string& job)
{
    SimDAG=&DAG;
    thresold=pthres;
    jobName=job;
    levels=new std::vector<std::vector<LandscapeNode*>*>();
    unusedGenes=std::list<int>();
    parentVector=std::vector<int>(SimDAG->size());
    for(int i=0; i<parentVector; i++){
        parentVector[i] = 0;
    }

    setupRoot();
    prepareGeneList();
    bool changed = true;
    if(changed && unusedGenes.size() > 0){
        changed = createNewLevel();
    }
}


LandscapeTree::setupRoot(){

   std::vector<LandscapeNode>* zeroLevel = new std::vector<LandscapeNode*>();
   LandscapeNode* root = new LandscapeNode();
   root->label = new std::string("Normal");
   root->familyNodesIDs = new std::vector<int>();
   root->nodeID = 0;
   root->prob = 1;
   root->parentID = -1;
   zeroLevel->push_back(root);
   levels->push_back(zeroLevel);

}

void LandscapeTree::prepareGeneList(std::vector<int>& unusedG ){

    for(int i=0; i<SimDAG->size() ; i++){   // each gene has an ID > 0 (0 is root's ID)
        unusedG.push_back(i+1);
    }

}

bool LandscapeTree::createNewLevel(){
    std::list<int> justUsedGenes = list<int>();
    std::vector<LandscapeNode*>* newLevel = new std::vector<LandscapeNode*>();
    for( LandscapeNode* node: levels[levels->size()-1]){
        for( int i=0 ; i<justUsedGenes.size(); i++){
            int gene = justUsedGenes[i];
            LandscapeNode* newNode;
            bool isInserted = computeNodeProbability(gene, node, newNode);  //newNode is passed by reference
            if(isInserted){
                newLevel->push_back(newNode);
            }
        }
        for( int i=0 ; i<unusedGenes.size(); i++){
            int gene = unusedGenes[i];
            LandscapeNode* newNode;
            bool isInserted = computeNodeProbability(gene, node, newNode);  //newNode is passed by reference
            if(isInserted){
                newLevel->push_back(newNode);
                justUsedGenes.push_back(gene);
                unusedGenes.remove(i);
            }
        }
    }
    if(newLevel->size() == 0){
        delete(newLevel);
        return false;
    }

    levels->push_back(newLevel);
    return true;
}


bool LandscapeTree::computeNodeProbability(int gene, LandscapeNode* source, LandscapeNode* destination){
    /* populate parent vecor */
    for(int i: source->familyNodesIDs){
        parentVector[i]=1;
    }
    parentVector[destination->nodeID]=1;

    double prob = 1;
    configParam((*SimDAG), parentVector, prob);

    /* restore parent vecor */
    for(int i: source->familyNodesIDs){
        parentVector[i]=0;
    }
    parentVector[destination->nodeID]=0;

    if(prob){

    }
}

