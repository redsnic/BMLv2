/*
 *  Algorithm.cpp
 *  
 *
 *  Created by Navodit Misra.
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
 */

/*	This is the main algorithm that controls and delegates specific tasks to other programs.  */

#include "Algorithm.h"

using namespace std;

/**
 * @brief execute
 * Main BML execution procedure
 *
 * @param input           struct containing the following fields
 * job                    job ID for output's names
 * inputFile              input data matrix
 * numberOfThreads        number of threads to be used by the program
 * nTrees                 number of tree reseeds
 * pthres                 thresold for inferring paths
 * performBootstrap       true if bootsrap should be preformed (false otherwise)
 * nRep                   number of bootstrap replicates to be generated
 * autoCutoff             if cutoff should be automaticalli computed
 * cutoff                 cutoff value (forced)
 *
 * @return                status of the execution (0 OK, else error)
 */
int execute(struct Input* input){

    try {

        srand((long) time(nullptr));

		vector< vector<bool> > data;
		vector<SimVar> simDAG;

		/*Learn the Bayesan Network */

        if(!structLearn(data, simDAG, input->inputFile, input->nTrees, input->pthres, input->job, input->numberOfThreads, input->cutoff, input->autoCutoff)){
            delete(input);
            return -1;
        }

		/* Perform bootstrap */

        if (input->performBootstrap) {
            bootstrapAnalysis(data, simDAG, input->nRep, input->job, input->numberOfThreads);
		}
        delete(input);
		return 0;

	} catch (int ex) {
		exceptionHandler(ex);  //	Exception handler (defined in IOHandelr)
	}
    delete(input);
    return -1;
}
