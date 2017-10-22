#ifndef LANDSCAPETREE_H
#include "IOHelper.h"
#ifndef DS
#include "DataStructures.h"
#define DS
#endif
#define LANDSCAPETREE_H

struct LandscapeNode{
    std::string* label;
    int nodeID;
    std::vector<int>* familyNodesIDs;
    int parentID; // -1 if root
    float prob;
};

class LandscapeTree
{
private:
    std::vector<std::vector<LandscapeNode*>*>* levels;
    std::list<int> unusedGenes;
    vector<SimVar>* SimDAG;
    std::string jobName;
    std::vector<int> parentVector;
    double thresold;
    void setupRoot();
    bool createNewLevel();
    bool computeNodeProbability(int, LandscapeNode*, LandscapeNode*);
    void prepareGeneList();
public:
    LandscapeTree(std::vector<SimVar> &SimDAG, double pthres, string &job);
    writeOnDisk();
};

#endif // LANDSCAPETREE_H
