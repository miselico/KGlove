/*
 * RDF2Co_occurence.cpp
 *
 *  Created on: Nov 29, 2016
 *
 * This implementation is more performant, but rather convoluted. Use KGloVe.cpp instead.
 *
 *      Author: cochez
 */

#include "RDF2Co_occurenceComplex.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/flyweight/flyweight.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>                   // for std::pair
#include <vector>

#include "BCA.h"
#include "graph/LabeledGraph.h"
#include "MyMaxPriorityQueue.h"
#include "nTripleParser.h"
#include "utils.h"


typedef boost::flyweight<std::string> flyString;

using namespace std;

namespace BCAOrder {

class my_bitset: public boost::dynamic_bitset<> {

private:
	boost::dynamic_bitset<>::size_type lastFound;
	boost::dynamic_bitset<>::size_type lastSetToTrue;

public:

	my_bitset() :
			lastFound(0), lastSetToTrue(0) {

	}

	boost::dynamic_bitset<>::size_type find_any() {
		if ((*this)[lastSetToTrue]) {
			return lastSetToTrue;
		}
		this->lastFound = this->find_next(lastFound);
		if (lastFound != this->npos) {
			return lastFound;
		} else {
			lastFound = find_first();
			return lastFound;
		}

	}

	dynamic_bitset& setTrueAndRecord(size_type n) {
		this->set(n, true);
		this->lastSetToTrue = n;
		return *this;
	}

};

//does the vector contain all numbers from 0 till all.size()-1 ?
bool allNumbersIn(std::vector<unsigned int> all) {
	std::sort(all.begin(), all.end());
	for (unsigned int i = 0; i < all.size(); i++) {
		if (all[i] != i) {
			return false;
		}
	}
	return true;
}

class Node {
public:
	//TODO this also works with std::set instead of unordered. Check what is faster for large graphs.
	int ID;
	std::unordered_set<Node*> inedgeSources;
	std::unordered_set<Node*> outedgeDestination;

	Node(int ID) :
			ID(ID), inedgeSources(2), outedgeDestination(2) {			//Intentionally empty
	}

	int getInDeg() const {
		return inedgeSources.size();
	}
	int getOutDeg() const {
		return outedgeDestination.size();
	}

};

// directed graph (single directed edge between an ordered pair of nodes)
//This implementation is made specifically for determineBCVcomputeOrder
class MyGraph: public std::vector<Node> {
public:
	MyGraph(std::shared_ptr<const QuickGraph::LabeledGraph> baseGraph) {
		cout << currentTime() << "copying graph" << endl;
		const std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
		const int totalNodes = nodes.size();
		this->reserve(totalNodes);
		for (int i = 0; i < totalNodes; i++) {
			this->emplace_back(i);
		}
		cout << currentTime() << "copying graph - nodes copied" << endl;

		//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
		//We also remove all self edges
		//for each node
		for (auto originalNode = nodes.begin(); originalNode < nodes.end(); originalNode++) {
			int src = originalNode - nodes.begin();
			//for each out edge of that node
			for (auto originalEdge = originalNode->edges.begin(); originalEdge < originalNode->edges.end(); originalEdge++) {
				int dst = originalEdge->targetIndex;
				if (src != dst) {
					std::unordered_set<Node*>::value_type sourceNode = &(*this)[src];
					std::unordered_set<Node*>::value_type destinationNode = &(*this)[dst];

					//std::pair<boost::unordered::iterator_detail::c_iterator<boost::unordered::detail::ptr_node<Node*> >, bool> a;
					bool inserted;
					std::tie(std::ignore, inserted) = sourceNode->outedgeDestination.insert(destinationNode);

					if (inserted) {
						destinationNode->inedgeSources.insert(sourceNode);
					}
				}
			}
		}

		cout << currentTime() << "graph copied" << endl;
	}
};

///**
// * Attempts to determine the order in which most BCA computations can be reused.
// * This is using the heuristic that
// *
// * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
// * 2. if no such node exists, then the node with highest in degree is selected
// *
// * the input graph MUST have nodes indexed 0 till Nodes()-1
// *
// */
vector<unsigned int> determineBCAcomputeOrder(std::shared_ptr<const QuickGraph::LabeledGraph> baseGraph) {

	MyGraph prunable_graph(baseGraph);
	const int unsigned startSize = baseGraph->nodes.size();
	baseGraph = NULL; //making sure we do not accidentally write to the original

	vector<unsigned int> finalOrder;
	finalOrder.reserve(startSize);

	my_bitset zeroOutDegrees;
	zeroOutDegrees.resize(startSize, false);

	boost::dynamic_bitset<> todo;
	todo.resize(startSize, true);

	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (std::vector<Node>::const_iterator node = prunable_graph.cbegin(); node != prunable_graph.cend(); node++) {
		if (node->getOutDeg() == 0) {
			zeroOutDegrees[node->ID] = true;
		}
	}

	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;

	unsigned long int n;
	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
		zeroOutDegrees[n] = false;
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		Node& toBeRemoved = prunable_graph[n];

		for (std::unordered_set<Node*>::const_iterator dependant = toBeRemoved.inedgeSources.cbegin(); dependant != toBeRemoved.inedgeSources.cend(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			std::unordered_set<Node*>::value_type toBeRemovedAdress = &toBeRemoved;
			(*dependant)->outedgeDestination.erase(toBeRemovedAdress);
		}
		//finalize, delete node
		//the only things remaining are some bookkeeping:
#ifndef NDEBUG
		//probably not necessary anyhow
		toBeRemoved.inedgeSources.clear();
#endif
		todo[n] = false;
		finalOrder.push_back(n);
		if (finalOrder.size() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.size() << "/" << startSize << " done" << endl;
		}
	}

	cout << currentTime() << "After first fast phase, " << finalOrder.size() << "/" << startSize << " nodes are done, starting iterative phase" << endl;

	//now the more general case including loops is handled
	MyMaxPriorityQueue<int> highestInDegree;
	//set-up the  highestInDegree PQ,

	highestInDegree.Size();

	for (std::vector<Node>::const_iterator node = prunable_graph.cbegin(); node < prunable_graph.cend(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		if (todo[node->ID]) {
			double priority = node->getInDeg() + 1;
			highestInDegree.Insert(node->ID, priority);
		}
	}

	//algo start

	while (finalOrder.size() < startSize) {
		//break at least one loop
		int k = -1;
		do {
			assert(highestInDegree.Size() > 0);
			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
			k = highestInDegree.PopMax();
		} while (!todo[k]);		//this check is needed because the priorities for nodes removed in the loop below are not removed from the queueu, only marked as done
		//k will be removed add all nodes which will get a zero out degree to set
		Node& kNode = prunable_graph[k];

		for (std::unordered_set<Node*>::const_iterator dependant = kNode.inedgeSources.cbegin(); dependant != kNode.inedgeSources.cend(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			(*dependant)->outedgeDestination.erase(&kNode);
		}

		//update the priorities of all nodes k is pointing to
		for (std::unordered_set<Node*>::const_iterator dest = kNode.outedgeDestination.cbegin(); dest != kNode.outedgeDestination.cend(); dest++) {
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that dest.id != k
			highestInDegree.SetPriority((*dest)->ID, (*dest)->getInDeg() - 1 + 1);
			//start removing already
			(*dest)->inedgeSources.erase(&kNode);
		}

		//finalize, the removal is already done in the adjecancy lists
#ifndef NDEBUG
		//probably not necessary anyhow
		kNode.inedgeSources.clear();
		kNode.outedgeDestination.clear();
#endif
		todo[k] = false;
		finalOrder.push_back(k);
		if (finalOrder.size() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.size() << "/" << startSize << " done" << endl;
		}

		if (finalOrder.size() == startSize) {
			break;
		}

		//remove as much as possible without having to break a loop
		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
			zeroOutDegrees[n] = false;
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			Node& toBeRemoved = prunable_graph[n];

			for (std::unordered_set<Node*>::const_iterator dependant = toBeRemoved.inedgeSources.cbegin(); dependant != toBeRemoved.inedgeSources.cend(); dependant++) {
				//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
				if ((*dependant)->getOutDeg() == 1) {
					zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
				}
				//start removing already
				(*dependant)->outedgeDestination.erase(&toBeRemoved);
			}
			//finalize, delete node
			//the only things remaining are some bookkeeping:
#ifndef NDEBUG
			//probably not necessary anyhow
			toBeRemoved.inedgeSources.clear();
#endif
			todo[n] = false;
			finalOrder.push_back(n);
			if (finalOrder.size() % infoFrequency == 0) {
				cout << currentTime() << finalOrder.size() << "/" << startSize << " done" << endl;
			}
		}

	}

	cout << currentTime() << "All done" << finalOrder.size() << "/" << startSize << " done" << endl;

	assert(startSize == finalOrder.size());
	assert(todo.find_first() == todo.npos);
	assert(allNumbersIn(finalOrder));
	return finalOrder;
}

void ComputeBCAOrder(const string & inputFilename, const string & outputFileName) {
	std::pair<std::shared_ptr<QuickGraph::LabeledGraph>, std::unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputFilename);
	std::shared_ptr<QuickGraph::LabeledGraph> graph = graphAndNodeIndex.first;
	cout << currentTime() << " Now computing BCV compute order" << endl;
	vector<unsigned int> order = determineBCAcomputeOrder(graph);
	cout << currentTime() << "end determining BCV compute order, writing to file" << endl;

	ofstream myfile(outputFileName);
	for (vector<unsigned int>::iterator iter = order.begin(); iter < order.end(); iter++) {
		myfile << int(*iter) << '\n';
	}
	myfile.flush();
	if (!myfile.good()) {
		throw "WTF";
	}
	myfile.close();
	cout << currentTime() << "done writing" << endl;
}

vector<unsigned int> readBCAOrder(const string & precomputedBCAOrderFile, unsigned int expectedCount) {
	vector<unsigned int> result;
	result.reserve(expectedCount);
	ifstream myfile(precomputedBCAOrderFile);
	unsigned int next;
	while (myfile >> next) {
		result.push_back(next);
	}
	myfile.close();

	if (result.size() != expectedCount) {
		throw "Size of the precomputed BCA order file is not as expected";
	}

	return result;
}

}

namespace co_occurence_computer {
/**
 * typedefs compatible with the input expected by glove
 */
typedef double real;

typedef struct cooccur_rec {
	int word1;
	int word2;
	real val;
	cooccur_rec(int word1, int word2, real val) :
			word1(word1), word2(word2), val(val) {

	}
} CREC;

/**
 * Get a table which assigns a unique index to each (predicate, object) pair. For PBCA
 *
 * The index will be strictly greater as zero since that is what glove uses to skip an entry.
 */
//static unordered_map<pair<flyString, flyString>, int, pairhash> createPairedWordIndexTable(shared_ptr<QuickGraph::LabeledGraph> graph) {
//
//	unordered_map<pair<flyString, flyString>, int, pairhash> table;
//
//	int counter = 1;
//	for (auto nodeI = graph->nodes.begin(); nodeI < graph->nodes.end(); nodeI++) {
//		for (auto edgeI = nodeI->edges.begin(); edgeI < nodeI->edges.end(); edgeI++) {
//			flyString predicate = edgeI->label;
//			flyString object = graph->nodes[edgeI->targetIndex].label;
//			std::pair<flyString, flyString> p(predicate, object);
//
//			if (table.find(p) == table.end()) {
//				std::pair<std::pair<flyString, flyString>, int> q(p, counter);
//				table.insert(q);
//				counter++;
//			}
//		}
//	}
//	return table;
//}
/*
 * The wordIndexTable is indexed from 0 as it should be, but we need IDs which start from 1. This function abstracts this away
 */

static unsigned int graphIDToGloveID(unsigned int graphID) {
	return graphID + 1;
}

/**
 * For each unique predicate in the graph add a unique ID it to the returned hash.
 *
 * The returned vector contains the strings in the same order as their IDs
 *
 * The used IDs range from graph.Nodes() till graph.Nodes+(numberOfUniquePredicates=returnValue.size()) exclusive
 */
static pair<vector<flyString>, std::unordered_map<flyString, unsigned int>> computePredicateIDs(shared_ptr<QuickGraph::LabeledGraph> graph) {
	std::unordered_map<flyString, unsigned int> preds;
	vector<flyString> labels;
	unsigned int currentID = graph->nodes.size();

	for (auto nodeI = graph->nodes.begin(); nodeI < graph->nodes.end(); nodeI++) {
		for (auto edgeI = nodeI->edges.begin(); edgeI < nodeI->edges.end(); edgeI++) {
			flyString edgeLabel = edgeI->label;
			if (preds.find(edgeLabel) == preds.end()) {
				preds[edgeLabel] = currentID;
				currentID++;
				labels.push_back(edgeLabel);
			}
		}
	}

	return pair<vector<flyString>, unordered_map<flyString, unsigned int> >(labels, preds);
}

bool isEntity(shared_ptr<QuickGraph::LabeledGraph> graph, unsigned int candidateIndex) {
	std::vector<QuickGraph::Edge> & edges = graph->nodes[candidateIndex].edges;
	for (std::vector<QuickGraph::Edge>::iterator edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
		flyString predicate = edgeI->label;
		if (predicate == RDF_TYPE) {
			QuickGraph::Node & targetNode = graph->nodes[edgeI->targetIndex];
			flyString object = targetNode.label;
			if (object == OWL_THING) {
				return true;
			}
		}
	}
	return false;
}

class Co_occurenceComputer {
private:
	//currently preceded by _ to prevent name clashes in not yet adapted methods
	vector<unsigned int> _order;
	std::shared_ptr<QuickGraph::LabeledGraph> _weightedGraph;
	//pair<vector<string>, unordered_map<string, int>> predLabelsAndIDs;
	unordered_map<flyString, unsigned int> _predGraphIDs;
	vector<flyString> _predicateLabels;

public:

	Co_occurenceComputer(const string & inputGraphFileName, weigher::GraphWeigher& weighingStrategy, const bool removeLiterals) {
		pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraph(inputGraphFileName, removeLiterals);
		_weightedGraph = graphAndNodeIndex.first;
		reWeigh(weighingStrategy);
		pair<vector<flyString>, unordered_map<flyString, unsigned int> > predLabelsAndIDs = computePredicateIDs(_weightedGraph);
		_predGraphIDs = predLabelsAndIDs.second;
		_predicateLabels = predLabelsAndIDs.first;

		_order = BCAOrder::determineBCAcomputeOrder(_weightedGraph);
	}

	Co_occurenceComputer(const string & inputGraphFileName, string precomputedBCAOrderFile, weigher::GraphWeigher& weighingStrategy) {
		pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputGraphFileName);
		_weightedGraph = graphAndNodeIndex.first;

		pair<vector<flyString>, unordered_map<flyString, unsigned int> > predLabelsAndIDs = computePredicateIDs(_weightedGraph);
		_predGraphIDs = predLabelsAndIDs.second;
		_predicateLabels = predLabelsAndIDs.first;

		_order = BCAOrder::readBCAOrder(precomputedBCAOrderFile, _weightedGraph->nodes.size());

		this->reWeigh(weighingStrategy);
	}

	void reWeigh(weigher::GraphWeigher& weighingStrategy) {
		weighingStrategy.weigh(_weightedGraph);
	}

	/**
	 *
	 * Compute the BCA score for each pair in the graph under the given weighing strategy.
	 * Additionally, adds a score for each edge as well.
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 *
	 *
	 */
	void computeFrequenciesIncludingEdges(double bca_alpha, double bca_eps, FILE * glove_input_file_out, bool normalize, bool onlyEntities) {

		const int infoFrequency = _order.size() > 1000000 ? 100000 : 10000;
		vector<std::shared_ptr<CompactBCV>> bcvCache;
		bcvCache.resize(this->_weightedGraph->nodes.size(), 0);

		int counter = 0;
		for (vector<unsigned int>::iterator iter = _order.begin(); iter < _order.end(); iter++) {
			const unsigned int focusWordGraphID = *iter;
			//		//only take specific one:
			//		if (candidateNode.second != "<http://dbpedia.org/ontology/Province>"){
			//			continue;
			//		}

			if (onlyEntities) {
				if (!isEntity(this->_weightedGraph, focusWordGraphID)) {
					continue;
				}
			}
			shared_ptr<CompactBCV> combinedbcv = computeBCAIncludingEdgesCached(_weightedGraph, focusWordGraphID, bca_alpha, bca_eps, _predGraphIDs, bcvCache);
			const int focusWordGloveID = graphIDToGloveID(focusWordGraphID);

			if (normalize) {
				//Removes value for focusword + normalizes and write that value out immediately.
				double totalSum = 0.0;
				for (vector<pair<unsigned int, double>>::const_iterator iter = combinedbcv->values.cbegin(); iter != combinedbcv->values.cend(); iter++) {
					if (iter->first != focusWordGraphID) {
						totalSum += iter->second;
					}
				}
				for (vector<pair<unsigned int, double>>::const_iterator iter = combinedbcv->values.cbegin(); iter != combinedbcv->values.cend(); iter++) {
					if (iter->first != focusWordGraphID) {
						int contextWordGraphID = iter->first;
						int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
						double freq = iter->second;
						CREC crec = CREC(focusWordGloveID, contextWordGloveID, freq / totalSum);
						//CREC crec = CREC { word1:contextWordGloveID, word2:focusWordGloveID, val: freq };
						fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
					}
				}
				//combinedbcv.removeEntry(focusWordGraphID);
				//combinedbcv.normalizeInPlace();
			} else {
				for (vector<pair<unsigned int, double>>::const_iterator iter = combinedbcv->values.cbegin(); iter != combinedbcv->values.cend(); iter++) {
					int contextWordGraphID = iter->first;
					int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
					double freq = iter->second;
					CREC crec = CREC(focusWordGloveID, contextWordGloveID, freq);
					//CREC crec = CREC { word1:contextWordGloveID, word2:focusWordGloveID, val: freq };
					fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
				}
			}
			counter++;
			if ((counter % infoFrequency) == 0) {
				cout << currentTime() << "Processed " << counter << "/" << _order.size() << " BCV computations" << endl;
			}
		}

	}

	void writeVocabFileIncludingEdges(FILE * glove_vocab_file_out) {
		//still need to write all node labels to the vocab file

		std::vector<QuickGraph::Node> & nodes = this->_weightedGraph->nodes;
		for (std::vector<QuickGraph::Node>::iterator iter = nodes.begin(); iter != nodes.end(); iter++) {
			const string & label = iter->label.get();
			fprintf(glove_vocab_file_out, "%s nofr\n", label.c_str());
		}

		//still need to write all predicates to the vocab file

		for (vector<flyString>::iterator it = _predicateLabels.begin(); it < _predicateLabels.end(); it++) {
			fprintf(glove_vocab_file_out, "%s nofr\n", it->get().c_str());
		}
	}

	/**
	 *
	 * Compute the BCA score for each noe pair in the graph under the given weighing strategy. Ignores predicates.
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 */
	void computeFrequencies(double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out, bool normalize, bool onlyEntities) {
		cerr << "computeFrequencies is not yet tested thoroughly, check results with care";
		const int infoFrequency = _order.size() > 1000000 ? 100000 : 10000;
		unordered_map<unsigned int, BCV> bcvCache;
		int counter = 0;

		for (vector<unsigned int>::iterator iter = _order.begin(); iter < _order.end(); iter++) {
			const unsigned int focusWordGraphID = *iter;

			if (onlyEntities) {
				if (!isEntity(this->_weightedGraph, focusWordGraphID)) {
					continue;
				}
			}
			BCV bcv = computeBCACached(_weightedGraph, focusWordGraphID, bca_alpha, bca_eps, bcvCache);

			int focusWordGloveID = graphIDToGloveID(focusWordGraphID);
			if (normalize) {
				bcv.removeEntry(focusWordGraphID);
				bcv.normalizeInPlace();
			}

			for (unordered_map<unsigned int, double>::iterator iter = bcv.begin(); iter != bcv.end(); iter++) {
				int contextWordGraphID = iter->first;
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter->second;
				CREC crec = CREC(focusWordGloveID, contextWordGloveID, freq);
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}
			counter++;
			if ((counter % infoFrequency) == 0) {
				cout << currentTime() << "Processed " << counter << "/" << _order.size() << " BCV computations" << endl;
			}
		}
		//still need to write all node labels to the vocab file

		std::vector<QuickGraph::Node> & nodes = this->_weightedGraph->nodes;
		for (std::vector<QuickGraph::Node>::iterator iter = nodes.begin(); iter != nodes.end(); iter++) {
			fprintf(glove_vocab_file_out, "%s nofr\n", iter->label.get().c_str());
		}
	}

};

class Co_occurenceComputer_Ultimate {
private:
	//currently preceded by _ to prevent name clashes in not yet adapted methods
	vector<unsigned int> _order;
	std::shared_ptr<QuickGraph::LabeledGraph> _weightedGraph;
	//pair<vector<string>, unordered_map<string, int>> predLabelsAndIDs;
	unordered_map<flyString, unsigned int> _predGraphIDs;
	vector<flyString> _predicateLabels;
	std::shared_ptr<QuickGraph::LabeledGraph> _weightedReverseGraph;
	vector<unsigned int> _orderReverse;

public:

	Co_occurenceComputer_Ultimate(const string & inputGraphFileName, weigher::GraphWeigher& weighingStrategy, weigher::GraphWeigher & reverseWeighingStrategy, const bool removeLiterals) {
		pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraph(inputGraphFileName, removeLiterals);
		_weightedGraph = graphAndNodeIndex.first;

		pair<vector<flyString>, unordered_map<flyString, unsigned int>> predLabelsAndIDs = computePredicateIDs(_weightedGraph);
		_predGraphIDs = predLabelsAndIDs.second;
		_predicateLabels = predLabelsAndIDs.first;

		_order = BCAOrder::determineBCAcomputeOrder(_weightedGraph);

		//reverse graph
		_weightedReverseGraph = QuickGraph::reverseGraph(_weightedGraph);
		_orderReverse = BCAOrder::determineBCAcomputeOrder(_weightedReverseGraph);
		this->reWeigh(weighingStrategy, reverseWeighingStrategy);
	}

	/**
	 * Carfull using this constructor!! The graph will not be weighted. reweigh has to be called before using the graph.
	 */
	Co_occurenceComputer_Ultimate(const string & inputGraphFileName, const bool removeLiterals) {
		pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraph(inputGraphFileName, removeLiterals);
		_weightedGraph = graphAndNodeIndex.first;

		pair<vector<flyString>, unordered_map<flyString, unsigned int>> predLabelsAndIDs = computePredicateIDs(_weightedGraph);
		_predGraphIDs = predLabelsAndIDs.second;
		_predicateLabels = predLabelsAndIDs.first;

		_order = BCAOrder::determineBCAcomputeOrder(_weightedGraph);

		//reverse graph
		_weightedReverseGraph = QuickGraph::reverseGraph(_weightedGraph);
		_orderReverse = BCAOrder::determineBCAcomputeOrder(_weightedReverseGraph);
	}

	void reWeigh(weigher::GraphWeigher& weighingStrategy, weigher::GraphWeigher & reverseWeighingStrategy) {
		weighingStrategy.weigh(_weightedGraph);
		reverseWeighingStrategy.weigh(this->_weightedReverseGraph);
	}

	/**
	 *
	 * Compute the BCA score for each pair in the graph under the given weighing strategy.
	 * Additionally, adds a score for each edge as well.
	 *
	 * Furthermore, it also performs a reverse walk and adds the result of that to the BCVs
	 * The reverse walk can be performed according to a different weighing strategy
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 */
	void computeFrequenciesIncludingEdgesTheUltimate(double bca_alpha, double bca_eps, FILE * glove_input_file_out, bool normalize, bool onlyEntities) {
		const int infoFrequency = _order.size() > 1000000 ? 100000 : 10000;
		//unordered_map<int, BCV> bcvForwardCache;
		//TODO : check whether this size is correct
		vector<std::shared_ptr<CompactBCV>> bcvForwardCache(_order.size());

		{			//filling the forward cache
			int counter = 0;
			for (vector<unsigned int>::iterator iter = _order.begin(); iter < _order.end(); iter++) {
				const unsigned int focusWordGraphID = *iter;

				if (onlyEntities) {
					if (!isEntity(this->_weightedGraph, focusWordGraphID)) {
						continue;
					}
				}
				//we only want the side effect of the BCV being added to the cache!
				computeBCAIncludingEdgesCached(_weightedGraph, focusWordGraphID, bca_alpha, bca_eps, _predGraphIDs, bcvForwardCache);
				counter++;
				if ((counter % infoFrequency) == 0) {
					cout << currentTime() << "Processed " << counter << "/" << _order.size() << " BCV FORWARD computations " << endl;
				}
			}
		}
		{			//filling the backward cache and writing comined outcomes to disk. This simultaniously empties the forwardCache to save memory.
			vector<std::shared_ptr<CompactBCV>> bcvBackwardCache(_order.size());
			int backwardCounter = 0;
			for (vector<unsigned int>::iterator iter = _orderReverse.begin(); iter < _orderReverse.end(); iter++) {
				const unsigned int focusWordGraphID = *iter;

				if (onlyEntities) {
					//This needs to be the original graph! The reversed graph does not have the information about what are entities!
					if (!isEntity(this->_weightedGraph, focusWordGraphID)) {
						continue;
					}
				}
				shared_ptr<CompactBCV> backWardBCV = computeBCAIncludingEdgesCached(_weightedReverseGraph, focusWordGraphID, bca_alpha, bca_eps, _predGraphIDs, bcvBackwardCache);
				backwardCounter++;
				if ((backwardCounter % infoFrequency) == 0) {
					cout << currentTime() << "Processed " << backwardCounter << "/" << _order.size() << " BCV BACKWARD computations" << endl;
				}
				//combine forward and backward:
				std::shared_ptr<CompactBCV> compactForwardBCV = bcvForwardCache[focusWordGraphID];
				//remove from forward cache. This is not 100% necessary, but helps ensuring program correctness, and saves a bit of memory. Any forward entry must only be needed once.
				bcvForwardCache[focusWordGraphID] = 0;

				BCV forwardBCV(*compactForwardBCV);
				forwardBCV.add(*backWardBCV);

				BCV &combinedBCV = forwardBCV;

				if (normalize) {
					combinedBCV.removeEntry(focusWordGraphID);
					combinedBCV.normalizeInPlace();
				}
				const int focusWordGloveID = graphIDToGloveID(focusWordGraphID);

				for (std::unordered_map<unsigned int, double>::const_iterator iter = combinedBCV.cbegin(); iter != combinedBCV.cend(); iter++) {
					int contextWordGraphID = iter->first;
					int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
					double freq = iter->second;
					CREC crec = CREC(focusWordGloveID, contextWordGloveID, freq);
					fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
				}
			}
		}
	}

	void writeVocabFileIncludingEdges(FILE * glove_vocab_file_out) {
		//still need to write all node labels to the vocab file

		std::vector<QuickGraph::Node> & nodes = this->_weightedGraph->nodes;
		for (std::vector<QuickGraph::Node>::iterator iter = nodes.begin(); iter != nodes.end(); iter++) {
			const string & label = iter->label.get();
			fprintf(glove_vocab_file_out, "%s nofr\n", label.c_str());
		}

		//still need to write all predicates to the vocab file

		for (vector<flyString>::iterator it = _predicateLabels.begin(); it < _predicateLabels.end(); it++) {
			fprintf(glove_vocab_file_out, "%s nofr\n", it->get().c_str());
		}
	}

};

//
///////////////////////////not yet adapted (reusing parts + optiimzations) static methods////////////////////
//
//
///**
// * Implementation incomplete!!!
// */
//void computeFrequenciesPushBCA(string filename, weigher::GraphWeigher& weighingStrategy, FILE *fout) {
//	cerr << "computeFrequenciesPushBCA is not yet completely implemented";
//	throw "not implemented";
//
//	pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, int> > graphAndNodeIndex = n3parser::buildRDFGraph(filename);
//	shared_ptr<QuickGraph::LabeledGraph> graph = graphAndNodeIndex.first;
//
//	unordered_map<string, int> wordIndexTable = graphAndNodeIndex.second;
//
//	unordered_map<pair<string, string>, int> pairwordIndexTable = createPairedWordIndexTable(graph);
//	std::shared_ptr<QuickGraph::LabeledGraph> weightedGraph = weighingStrategy.weigh(graph);
//
//	graph.Clr();
//
//	cout << "done weighing" << endl;
//
//	cerr << "TODO : check - does the indexing for glove have to start from 1 or 0 ??";
//	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {
//		if (weightedGraph->GetNDat(i).SearchCh('<') != 0) {
//			continue;
//		}
//		PBCV bcv = computePBCA(weightedGraph, i, 0.10, 0.000001);
//		int subjectIndex = wordIndexTable.GetDat(weightedGraph->GetNDat(i));
//		for (unordered_map<pair<int, int>, double>::iterator iter = bcv.begin(); iter < bcv.end(); iter++) {
//			WeightedPredicate wpred = weightedGraph->GetEDat(iter.first.first);
//			string pred = wpred.P();
//			string obj = weightedGraph->GetNDat(iter.first.second);
//			int offset = wordIndexTable.size();
//			int pairIndex = pairwordIndexTable.GetDat(pair<string, string>(pred, obj)) + offset;
//			double freq = iter.second;
//
//			CREC crec = CREC { word1:subjectIndex, word2:pairIndex, val: freq };
//			fwrite(&crec, sizeof(CREC), 1, fout);
//
//			cout << subjectIndex << " has " << pairIndex << " freq " << freq << endl;
//		}
//	}
//
//	//	TTmStopWatch w (true);
//	//	int needed = 10000;
//	//	int * selected = (int*) malloc(needed * sizeof(int));
//	//	int skipped = 0;
//	//	for (int i = 0; i < needed + skipped && i < weightedGraph->GetNodes(); ++i) {
//	//		if (weightedGraph->GetNDat(i).SearchCh('<') != 0) {
//	//			++skipped;
//	//			continue;
//	//		} else {
//	//			selected[i - skipped] = i;
//	//		}
//	//	}
//	//	for (int index = 0; index < needed && index < weightedGraph->GetNodes(); index++) {
//	//		int id = selected[index];
//	//		PBCV bcv = computePBCA(weightedGraph, id, 0.10, 0.000000000001);
//	//		//cout << bcv.toString(weightedGraph) << endl;
//	//		if (index % 1000 == 0) {
//	//			cout << "another 1000" << weightedGraph->GetNDat(id).CStr() << "->" << bcv.toString(weightedGraph) << endl;
//	//		}
//	//	}
//	//
//	//	w.Stop();
//	//	cout << w.GetMSecInt() << "ms" << endl;
//	return;
//}
//
}//end namespace co_occurence_computer

namespace RDF2CO {

//precompute the BCA order
//void performExperiments() {
//	string graphInputFile = "../../datasets/dbPedia/allData27_30M.nt";
//	auto BCAOrderFile = "BCAComputeOrder.txt";
//	BCAOrder::ComputeBCAOrder(graphInputFile, BCAOrderFile);
//}

void performExperiments() {

	cerr << "TODO: check whether increasing locality (instead of followig chains) in BCA order improves performance" << endl;
	cerr << "note: corrected fixing of paint on edges. Multiplied with alpha before fixing!" << endl;
	cerr << "TODO: check the following: at each BCV computation a bit of page rank is dropped. When reused, each time this bit is multiplied. Does this mean that in the end a lot of pagerank is lost?"
			<< endl;
	cerr << "enable -NDEBUG for actual runs" << endl;
//string graphInputFile = "../../datasets/dbPedia/allData27_30M.nt";
//string graphInputFile = "../../datasets/wikidata-simple-statements-10_000000-sample.nt";

//string graphInputFile = "SmallTest.nt";
//string graphInputFile = "SmallTest9_loop.nt";
//string graphInputFile = "SmallTest11_simpleLoop.nt";
	string graphInputFile = "proteinInteraction.nt";
	weigher::UniformWeigher weigher;
//two options, if the BCA order has been precomputed:
//string precomputedBCAOrderFile = "BCAComputeOrder.txt";
//Co_occurenceComputer c(graphInputFile, precomputedBCAOrderFile, weigher);
//if it has not been precomputed:
	co_occurence_computer::Co_occurenceComputer c(graphInputFile, weigher, false);
//	co_occurence_computer::Co_occurenceComputer_Ultimate c(graphInputFile, weigher, weigher);
//now, c can be used to compute co_occurence matrices

	string glove_vocab_file = "glove_vocab_file_" + boost::replace_all_copy(graphInputFile, ".", "-") + ".txt";
	FILE* glove_vocab_file_out = fopen(glove_vocab_file.c_str(), "w");
	c.writeVocabFileIncludingEdges(glove_vocab_file_out);
	fclose(glove_vocab_file_out);
	cout << "vocab written to " << glove_vocab_file << endl;

	for (double alpha = 0.1; alpha <= 0.5; alpha += 0.1) {
		for (double eps = 0.00001; eps <= 0.0001; eps *= 10) {

			string glove_input_file = "glove_input_file_out_alpha_" + boost::lexical_cast<std::string>(alpha) + "_eps_" + boost::lexical_cast<std::string>(eps) + "_file_"
					+ boost::replace_all_copy(graphInputFile, ".", "-") + ".bin";

			cout << "writing to " << glove_input_file << endl;

			FILE* glove_input_file_out = fopen(glove_input_file.c_str(), "w");

//			//c.computeFrequenciesIncludingEdgesTheUltimate(alpha, eps, glove_input_file_out, glove_vocab_file_out, true, false);

			c.computeFrequenciesIncludingEdges(alpha, eps, glove_input_file_out, false, false);

			fclose(glove_input_file_out);
		}
	}
}


void ParameterizedRun::parametrizedUltimateRun(Parameters& param) {
	param.check();

	int outputFileNumber = 0;
	for (std::vector<std::tuple<std::string, bool, bool>>::const_iterator graphs = param.graphs.begin(); graphs != param.graphs.end(); graphs++) {
		const std::string graphInputFile = std::get<0>(*graphs);
		const bool removeLiterals = std::get<1>(*graphs);
		const bool saveVocabulary = std::get<2>(*graphs);

		co_occurence_computer::Co_occurenceComputer_Ultimate c(graphInputFile, removeLiterals);
		if (saveVocabulary) {
			string glove_vocab_file_out_name = "glove_vocab_file_" + boost::replace_all_copy(graphInputFile, ".", "-") + ".txt";
			FILE* glove_vocab_file_out = fopen(glove_vocab_file_out_name.c_str(), "w");
			if (!glove_vocab_file_out) {
				throw "vocab file could not be opened";
			}
			c.writeVocabFileIncludingEdges(glove_vocab_file_out);
			fclose(glove_vocab_file_out);
		}
		string fileNamePART1 = boost::replace_all_copy(graphInputFile, ".", "_") + "-" + (removeLiterals ? "no_literals" : "with_literals");

		for (std::vector<std::pair<weigher::GraphWeigher&, weigher::GraphWeigher&>>::const_iterator weighers = param.weighers.begin(); weighers != param.weighers.end(); weighers++) {
			weigher::GraphWeigher & forwardWeigher = weighers->first;
			weigher::GraphWeigher & reverseWeigher = weighers->second;
			c.reWeigh(forwardWeigher, reverseWeigher);
			string fileNamePART2 = fileNamePART1 + "-forward_" + forwardWeigher.getName() + "-reverse_" + reverseWeigher.getName();
			for (std::vector<double>::const_iterator alphas = param.alphas.cbegin(); alphas != param.alphas.cend(); alphas++) {
				double alpha = *alphas;
				for (std::vector<double>::const_iterator epss = param.epss.cbegin(); epss != param.epss.cend(); epss++) {
					double eps = *epss;
					for (std::vector<bool>::const_iterator normalizes = param.normalize.cbegin(); normalizes != param.normalize.cend(); normalizes++) {
						bool normalize = *normalizes;
						for (std::vector<bool>::const_iterator onlyEntitiess = param.onlyEntities.cbegin(); onlyEntitiess != param.onlyEntities.cend(); onlyEntitiess++) {
							bool onlyEntities = *onlyEntitiess;

							//creating filename
							string outputFile = "glove_input_file-" + fileNamePART2 + "-alpha_" + boost::lexical_cast<std::string>(alpha) + "-eps_" + boost::lexical_cast<std::string>(eps)
									+ (normalize ? "-normalize_yes" : "-normalize_no") + (onlyEntities ? "-onlyEntities_yes" : "-onlyEntities_no") + ".bin";
							FILE* glove_input_file_out = fopen(outputFile.c_str(), "w");
							if (!glove_input_file_out) {
								throw "coocurence file could not be written to";
							}
							c.computeFrequenciesIncludingEdgesTheUltimate(alpha, eps, glove_input_file_out, normalize, onlyEntities);
							fclose(glove_input_file_out);
							outputFileNumber++;
						}
					}
				}
			}
		}
	}
}

}			//end RDF2CO
