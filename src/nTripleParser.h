/*
 * nTripleParser.h
 * This is a rudimentary parser for n-triples. It does not properly check for errors on the input, so make sure to feed proper n-triples.
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 */

#ifndef NTRIPLEPARSER_H_
#define NTRIPLEPARSER_H_

#include <iostream>
#include <utility>
#include <unordered_map>

#include "graph/LabeledGraph.h"

namespace n3parser {

std::pair<std::shared_ptr<QuickGraph::LabeledGraph>, std::unordered_map<std::string, unsigned int> > buildRDFGraph(const std::string & filename);

std::pair<std::shared_ptr<QuickGraph::LabeledGraph>, std::unordered_map<std::string, unsigned int> > buildRDFGraphIgnoreLiterals(const std::string & filename);

std::pair<std::shared_ptr<QuickGraph::LabeledGraph>, std::unordered_map<std::string, unsigned int> > buildRDFGraph(const std::string & filename, const bool removeLiterals);


}

#endif /* NTRIPLEPARSER_H_ */
