/*
 * InitTree.h
 *
 *  Created on: 24 ago 2017
 *      Author: rossi
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>

#ifndef INITTREE_H_

#ifndef DS
#include "DataStructures.h"
#define DS
#endif
#include "Utils.h"

void sanityCheckTree(vector<NODE>& Tree);
void initTree( vector<NODE>& Tree, vector<FAM>& BN);

#define DEBUG false

#define INITTREE_H_





#endif /* INITTREE_H_ */
