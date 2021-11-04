/*
 * RDF2Co_occurence.h
 *
 *  Created on: Nov 29, 2016
 *
 *  Major overhaul May 2019: removed many options to run more complex versions of the co-occurence matrix creation. they are still possible using the coocurence tools.
 *  This is less performant, but uses a bit less memory and the code becomes more easy to understand and experiment with.
 *
 *      Author: cochez
 */

#ifndef KGLOVE_HA_
#define KGLOVE_HA_

#include <string>
#include "GraphWeigher.h"
#include <tuple>

namespace KGloVe {

class Parameters {

public:
	/**
	 * The name of the graph, removeLiterals, should the vocabulary file be stored
	 */
	std::vector<std::tuple<std::string, bool, bool>> graphs;
	/**
	 * Weighers to be used. Note, only forward weighing is supported in this implementation.
	 * To get reverse weighing, feed the graph as desired and merge the different outcomes using the co-occurrence tools.
	 */
	std::vector<weigher::GraphWeigher *> weighers;
	//Values for alpha (BCA)
	std::vector<double> alphas;
	//Values for epsilon (BCA)
	std::vector<double> epss;
	//Only create scores for entities in the co-occurence matrix?
	std::vector<bool> onlyEntities;
	//Should edges also get a score in the co-occurence matrix?
	std::vector<bool> includeEdges;

	std::string outputPrefix;

	void check() {
		for (auto graph : graphs) {
			const std::string graphInputFile = std::get<0>(graph);
			// This check was needed because paths were not dealth with properly in earlier versions. Should be ok now.
//			if (graphInputFile.find('/') != std::string::npos) {
//				throw "Input path cannot contain '/'.";
//			}
		}

		if (!(graphs.size() > 0)) {
			throw "no graphs";
		}
		if (!(weighers.size() > 0)) {
			throw "no weighers";
		}
		if (!(alphas.size() > 0)) {
			throw "no alphas";
		}
		if (!(epss.size() > 0)) {
			throw "no epss";
		}
		if (!(onlyEntities.size() == 1 || (onlyEntities.size() == 2 && onlyEntities[0] != onlyEntities[1]))) {
			throw "onlyEntities can only be [true], [false], [true, false], or [false,true]";
		}
		if (!(includeEdges.size() == 1 || (includeEdges.size() == 2 && includeEdges[0] != includeEdges[1]))) {
			throw "onlyEntities can only be [true], [false], [true, false], or [false,true]";
		}
	}

};

/**
 * Will perform an minimal run on each permutation of the parameters.
 * This saves the co-occurence matrices from a simple forward BCA computation for each paramter permutation.
 *
 */
void parametrizedRun(Parameters & param);

}
#endif /* RDF2CO_OCCURENCE_HA_ */
