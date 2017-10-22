/*
 *  InitTree.cpp
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

#include "InitTree.h"

/**
 * Checks if there are problems with the created tree, if there are errors an exception is thrown.
 * This is a debug procedure.
 *
 * @param Tree           the tree to be checked
 */
void sanityCheckTree(vector<NODE>& Tree) {

	if (DEBUG) {
		if (Tree[0].children.size() != 1) {
			cout << "PROBLEM in InitTree.cpp. Zero has "
					<< Tree[0].children.size() << " children";
			throw 2;
		}
		for (unsigned int i = 1; i < Tree.size(); i++) {
			unsigned int parent = Tree[i].parent;
			if ((Tree[parent].children[0] != i) && (Tree[parent].children[1] != i)) {
				cout << "PROBLEM in InitTree.cpp";
				throw 2;
			}
		}
	}

}

/**
 * Reconstructs the phylogenetic tree by adding the internal nodes and the edges.
 *
 * @param tree           A phylogenetic tree with only vertex initialized (modified).
 *                       The new internal nodes will be added (after the leaves)
 *                       and edges will be created by randomly grouping nodes
 *                       and setting the hypotetical genome of internal nodes
 *                       so that there will be a number of mutations lower than the
 *                       number of observed mutations.
 *                       (as a mutations cannot be reversed by hypothesis)
 *                       The number of mutations for each incoming edge of a node will be also computed and saved.
 * @param bayesanNet     The bayesan network (with globalPruning) (unmodified)
 */
void initTree(vector<NODE>& tree, vector<FAM>& bayesanNet) {

	vector<int> Orphans;
    int numberOfSamples = int(tree.size());

	for (int i = 1; i < numberOfSamples; i++) {                         // at first each node is considered as orphan
		Orphans.push_back(i);
	}

    vector<int> numberOfMutations(0);                                      // the number of mutations seen for each node is stored
	for (unsigned int i = 0; i < bayesanNet.size(); i++) {
		numberOfMutations.push_back(bayesanNet[i].numberOfMutations);   // (the data is collected by the bayesanNet)
	}

	/* creation of internal nodes */

	for (int i = numberOfSamples; i < 2 * numberOfSamples - 2; i++) {   // being a complete binary tree the leaves are n+1 if n are the internal nodes

        double y;
        do{
            y = (double) rand() / RAND_MAX;
        }while(y>=1);
        int x = int(y * Orphans.size());                                // takes two random orphans
		int rightChild = Orphans[x];

		Orphans.erase(Orphans.begin() + x);

        do{
            y = (double) rand() / RAND_MAX;
        }while(y>=1);

		x = int(y * Orphans.size());
		int leftChild = Orphans[x];

		Orphans.erase(Orphans.begin() + x);

		/* creates an internal node */

		NODE n;                                                          // sets them as left and right children of a new node
		n.children.push_back(rightChild);
		n.children.push_back(leftChild);

		for(unsigned int k = 0; k < tree[0].genotype.size(); k++)        // each time a new internal vertex is created, check each gene of the samples
		{
			bool gbit = tree[rightChild].genotype[k] * tree[leftChild].genotype[k];  // true if both are mutated

			if (numberOfMutations[k] > numberOfSamples - 2) {            // internal nodes can't introduce unobserved mutations
				gbit = false;
			}

			n.genotype.push_back(gbit);                                  // create the predicted genotype for the internal node

			if (gbit == true) {										     // increase the number of mutations of this new internal node
				numberOfMutations[k]++;
			}
		}

		n.numberOfDescendantSamples = tree[rightChild].numberOfDescendantSamples + tree[leftChild].numberOfDescendantSamples;
		tree.push_back(n);                                               // add the new node to the tree

		tree[rightChild].parent = i;
		tree[leftChild].parent = i;


		Orphans.push_back(tree.size() - 1);

	}

	if (Orphans.size() != 1) {  // there sould be a single root
		cout << "Problem initializing tree. " << Orphans.size() << '\n';
	} else {                    // set root at index 0
		int rightChild = tree.size() - 1;
		tree[0].children.push_back(rightChild);
		tree[rightChild].parent = 0;
	}

	/* computes the number of mutations for all incoming edges of each node */

	for (unsigned int i = 1; i < tree.size(); i++) {
		for (unsigned int j = 0; j < tree[i].genotype.size(); j++) {
			if (tree[i].genotype[j] != tree[tree[i].parent].genotype[j]) {
				tree[i].numberOfMutationsOnEdge.push_back(j);
			}
		}
	}

	sanityCheckTree(tree);

}
