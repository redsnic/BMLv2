#ifndef ALGORITHM_CPP
#define ALGORITHM_CPP

/* Generic includes */
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>

#ifndef DS
#include "DataStructures.h"
#define DS
#endif

#include "IOHelper.h"
#include "GenerateReplicate.h"
#include "Utils.h"
#include "InitTree.h"
#include "ParseDataMatrix.h"
#include "SearchDAGs.h"
#include "SearchTrees.h"
#include "StructLearn.h"
#include "ui_mainmenu.h"

/* Data Structures & Utils */

int execute(Input* input);

#endif // ALGORITHM_CPP



