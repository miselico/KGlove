/*
 * RDF2Co_occurence.cpp
 *
 *  Created on: Nov 29, 2016
 *
 *  Major overhaul May 2019: derived from the more complex RDF2Co_occurenceComples.cpp removed many options to run more complex versions of the co-occurence matrix creation. They are still possible using the coocurence matrix tools.
 *  This is less performant, but uses a bit less memory and the code becomes more easy to understand and experiment with.
 *
 *  This also removes the option to run pushBCA
 *      Author: cochez
 */

#include <boost/algorithm/string/replace.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/flyweight/flyweight.hpp>
#include <boost/endian.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>                   // for std::pair
#include <vector>

#include "BCA.h"
#include "graph/LabeledGraph.h"
#include "GraphWeigher.h"
#include "MyMaxPriorityQueue.h"
#include "nTripleParser.h"
#include "KGloVe.h"
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
			ID(ID), inedgeSources(2), outedgeDestination(2) { //Intentionally empty
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
		const std::vector<QuickGraph::Node> &nodes = baseGraph->nodes;
		const int totalNodes = nodes.size();
		this->reserve(totalNodes);
		for (int i = 0; i < totalNodes; i++) {
			this->emplace_back(i);
		}
		cout << currentTime() << "copying graph - nodes copied" << endl;

		//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
		//We also remove all self edges
		//for each node
		for (auto originalNode = nodes.begin(); originalNode < nodes.end();
				originalNode++) {
			int src = originalNode - nodes.begin();
			//for each out edge of that node
			for (auto originalEdge = originalNode->edges.begin();
					originalEdge < originalNode->edges.end(); originalEdge++) {
				int dst = originalEdge->targetIndex;
				if (src != dst) {
					std::unordered_set<Node*>::value_type sourceNode =
							&(*this)[src];
					std::unordered_set<Node*>::value_type destinationNode =
							&(*this)[dst];

					//std::pair<boost::unordered::iterator_detail::c_iterator<boost::unordered::detail::ptr_node<Node*> >, bool> a;
					bool inserted;
					std::tie(std::ignore, inserted) =
							sourceNode->outedgeDestination.insert(
									destinationNode);

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
vector<unsigned int> determineBCAcomputeOrder(
		std::shared_ptr<const QuickGraph::LabeledGraph> baseGraph) {

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
	for (std::vector<Node>::const_iterator node = prunable_graph.cbegin();
			node != prunable_graph.cend(); node++) {
		if (node->getOutDeg() == 0) {
			zeroOutDegrees[node->ID] = true;
		}
	}

	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;

	unsigned long int n;
	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
		zeroOutDegrees[n] = false;
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		Node &toBeRemoved = prunable_graph[n];

		for (std::unordered_set<Node*>::const_iterator dependant =
				toBeRemoved.inedgeSources.cbegin();
				dependant != toBeRemoved.inedgeSources.cend(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			std::unordered_set<Node*>::value_type toBeRemovedAdress =
					&toBeRemoved;
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
			cout << currentTime() << finalOrder.size() << "/" << startSize
					<< " done" << endl;
		}
	}

	cout << currentTime() << "After first fast phase, " << finalOrder.size()
			<< "/" << startSize << " nodes are done, starting iterative phase"
			<< endl;

	//now the more general case including loops is handled
	MyMaxPriorityQueue<int> highestInDegree;
	//set-up the  highestInDegree PQ,

	highestInDegree.Size();

	for (std::vector<Node>::const_iterator node = prunable_graph.cbegin();
			node < prunable_graph.cend(); node++) {
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
		} while (!todo[k]);	//this check is needed because the priorities for nodes removed in the loop below are not removed from the queueu, only marked as done
		//k will be removed add all nodes which will get a zero out degree to set
		Node &kNode = prunable_graph[k];

		for (std::unordered_set<Node*>::const_iterator dependant =
				kNode.inedgeSources.cbegin();
				dependant != kNode.inedgeSources.cend(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			(*dependant)->outedgeDestination.erase(&kNode);
		}

		//update the priorities of all nodes k is pointing to
		for (std::unordered_set<Node*>::const_iterator dest =
				kNode.outedgeDestination.cbegin();
				dest != kNode.outedgeDestination.cend(); dest++) {
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that dest.id != k
			highestInDegree.SetPriority((*dest)->ID,
					(*dest)->getInDeg() - 1 + 1);
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
			cout << currentTime() << finalOrder.size() << "/" << startSize
					<< " done" << endl;
		}

		if (finalOrder.size() == startSize) {
			break;
		}

		//remove as much as possible without having to break a loop
		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
			zeroOutDegrees[n] = false;
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			Node &toBeRemoved = prunable_graph[n];

			for (std::unordered_set<Node*>::const_iterator dependant =
					toBeRemoved.inedgeSources.cbegin();
					dependant != toBeRemoved.inedgeSources.cend();
					dependant++) {
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
				cout << currentTime() << finalOrder.size() << "/" << startSize
						<< " done" << endl;
			}
		}

	}

	cout << currentTime() << "All done" << finalOrder.size() << "/" << startSize
			<< " done" << endl;

	assert(startSize == finalOrder.size());
	assert(todo.find_first() == todo.npos);
	assert(allNumbersIn(finalOrder));
	return finalOrder;
}

void ComputeBCAOrder(const string &inputFilename,
		const string &outputFileName) {
	std::pair<std::shared_ptr<QuickGraph::LabeledGraph>,
			std::unordered_map<string, unsigned int> > graphAndNodeIndex =
			n3parser::buildRDFGraphIgnoreLiterals(inputFilename);
	std::shared_ptr<QuickGraph::LabeledGraph> graph = graphAndNodeIndex.first;
	cout << currentTime() << " Now computing BCV compute order" << endl;
	vector<unsigned int> order = determineBCAcomputeOrder(graph);
	cout << currentTime()
			<< "end determining BCV compute order, writing to file" << endl;

	ofstream myfile(outputFileName);
	for (vector<unsigned int>::iterator iter = order.begin();
			iter < order.end(); iter++) {
		myfile << int(*iter) << '\n';
	}
	myfile.flush();
	if (!myfile.good()) {
		throw "WTF";
	}
	myfile.close();
	cout << currentTime() << "done writing" << endl;
}

vector<unsigned int> readBCAOrder(const string &precomputedBCAOrderFile,
		unsigned int expectedCount) {
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
typedef int64_t real_compatible_int;

typedef struct cooccur_rec {
	int word1;
	int word2;
	real val;
	cooccur_rec(int word1, int word2, real val) :
			word1(word1), word2(word2), val(val) {

	}
} CREC;

void writeCrec(CREC *crec, FILE *toFile) {
	int word1 = boost::endian::native_to_little(crec->word1);
	int word2 = boost::endian::native_to_little(crec->word2);

	fwrite(&word1, sizeof(int), 1, toFile);
	fwrite(&word2, sizeof(int), 1, toFile);

	static_assert(sizeof(real) == sizeof(real_compatible_int),"Endian conversion requires the same number of bytes in real and real_compatible_int");
	real float_val = crec->val;

	// from https://stackoverflow.com/a/20452215
	const int __one__ = 1;
	const bool isCpuLittleEndian = 1 == *(char*) (&__one__); // CPU endianness
	const bool isFileLittleEndian = true;  // output endianness

	if (isCpuLittleEndian ^ isFileLittleEndian) {
		char *pDouble = (char*) (double*) (&float_val);
		for (unsigned int i = 0; i < sizeof(real_compatible_int); ++i) {
			char byte = pDouble[sizeof(real_compatible_int) - 1 - i];
			fwrite(&byte, sizeof(char), 1, toFile);
		}
	} else {
		fwrite(&float_val, sizeof(real), 1, toFile);
	}
}

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
static pair<vector<flyString>, std::unordered_map<flyString, unsigned int>> computePredicateIDs(
		shared_ptr<QuickGraph::LabeledGraph> graph) {
	std::unordered_map<flyString, unsigned int> preds;
	vector<flyString> labels;
	unsigned int currentID = graph->nodes.size();

	for (auto nodeI = graph->nodes.begin(); nodeI < graph->nodes.end();
			nodeI++) {
		for (auto edgeI = nodeI->edges.begin(); edgeI < nodeI->edges.end();
				edgeI++) {
			flyString edgeLabel = edgeI->label;
			if (preds.find(edgeLabel) == preds.end()) {
				preds[edgeLabel] = currentID;
				currentID++;
				labels.push_back(edgeLabel);
			}
		}
	}

	return pair<vector<flyString>, unordered_map<flyString, unsigned int> >(
			labels, preds);
}

bool isEntity(shared_ptr<QuickGraph::LabeledGraph> graph,
		unsigned int candidateIndex) {
	std::vector<QuickGraph::Edge> &edges = graph->nodes[candidateIndex].edges;
	for (std::vector<QuickGraph::Edge>::iterator edgeI = edges.begin();
			edgeI != edges.end(); edgeI++) {
		flyString predicate = edgeI->label;
		if (predicate == RDF_TYPE) {
			QuickGraph::Node &targetNode = graph->nodes[edgeI->targetIndex];
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
	unordered_map<flyString, unsigned int> _predGraphIDs;
	vector<flyString> _predicateLabels;

public:

	Co_occurenceComputer(const string &inputGraphFileName,
			weigher::GraphWeigher &weighingStrategy,
			const bool removeLiterals) {
		pair<shared_ptr<QuickGraph::LabeledGraph>,
				unordered_map<string, unsigned int> > graphAndNodeIndex =
				n3parser::buildRDFGraph(inputGraphFileName, removeLiterals);
		_weightedGraph = graphAndNodeIndex.first;
		reWeigh(weighingStrategy);

		pair<vector<flyString>, unordered_map<flyString, unsigned int> > predLabelsAndIDs =
				computePredicateIDs(_weightedGraph);
		_predGraphIDs = predLabelsAndIDs.second;
		_predicateLabels = predLabelsAndIDs.first;

		_order = BCAOrder::determineBCAcomputeOrder(_weightedGraph);
	}

	/**
	 * Warning: no weighing startegy is provides, so no weiging has happened after this constructor reWeigh(weigher::GraphWeigher) must be called
	 */
	Co_occurenceComputer(const string &inputGraphFileName,
			const bool removeLiterals) {
		pair<shared_ptr<QuickGraph::LabeledGraph>,
				unordered_map<string, unsigned int> > graphAndNodeIndex =
				n3parser::buildRDFGraph(inputGraphFileName, removeLiterals);
		_weightedGraph = graphAndNodeIndex.first;

		pair<vector<flyString>, unordered_map<flyString, unsigned int> > predLabelsAndIDs =
				computePredicateIDs(_weightedGraph);
		_predGraphIDs = predLabelsAndIDs.second;
		_predicateLabels = predLabelsAndIDs.first;

		_order = BCAOrder::determineBCAcomputeOrder(_weightedGraph);
	}

//	Co_occurenceComputer(const string & inputGraphFileName, string precomputedBCAOrderFile, weigher::GraphWeigher& weighingStrategy) {
//		pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputGraphFileName);
//		_weightedGraph = graphAndNodeIndex.first;
//
//		pair<vector<flyString>, unordered_map<flyString, unsigned int> > predLabelsAndIDs = computePredicateIDs(_weightedGraph);
//		_predGraphIDs = predLabelsAndIDs.second;
//		_predicateLabels = predLabelsAndIDs.first;
//
//		_order = BCAOrder::readBCAOrder(precomputedBCAOrderFile, _weightedGraph->nodes.size());
//
//		this->reWeigh(weighingStrategy);
//	}

	void reWeigh(weigher::GraphWeigher &weighingStrategy) {
		weighingStrategy.weigh(_weightedGraph);
	}

	/**
	 *
	 * Compute the BCA score for each node pair in the graph under the given weighing strategy. Ignores predicates in the sense that they are not given a value in the co_occurence matrix
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 */
	void computeFrequenciesNoEdges(double bca_alpha, double bca_eps,
			FILE *glove_input_file_out, bool onlyEntities) {
		cerr
				<< "computeFrequencies is not yet tested thoroughly, check results with care";
		const int infoFrequency = _order.size() > 1000000 ? 100000 : 10000;
		unordered_map<unsigned int, BCV> bcvCache;
		int counter = 0;

		for (vector<unsigned int>::iterator iter = _order.begin();
				iter < _order.end(); iter++) {
			const unsigned int focusWordGraphID = *iter;

			BCV bcv = computeBCACached(_weightedGraph, focusWordGraphID,
					bca_alpha, bca_eps, bcvCache);

			if (onlyEntities) {
				if (!isEntity(this->_weightedGraph, focusWordGraphID)) {
					continue;
				}
			}

			int focusWordGloveID = graphIDToGloveID(focusWordGraphID);

			for (unordered_map<unsigned int, double>::iterator iter =
					bcv.begin(); iter != bcv.end(); iter++) {
				int contextWordGraphID = iter->first;
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter->second;
				CREC crec = CREC(focusWordGloveID, contextWordGloveID, freq);
				writeCrec(&crec, glove_input_file_out);
			}
			counter++;
			if ((counter % infoFrequency) == 0) {
				cout << currentTime() << "Processed " << counter << "/"
						<< _order.size() << " BCV computations" << endl;
			}
		}

	}

	void writeVocabFileNoEdges(FILE *glove_vocab_file_out) {
		//write all node labels to the vocab file

		std::vector<QuickGraph::Node> &nodes = this->_weightedGraph->nodes;
		for (std::vector<QuickGraph::Node>::iterator iter = nodes.begin();
				iter != nodes.end(); iter++) {
			const string &label = iter->label.get();
			fprintf(glove_vocab_file_out, "%s nofr\n", label.c_str());
		}
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
	void computeFrequenciesIncludingEdges(double bca_alpha, double bca_eps,
			FILE *glove_input_file_out, bool onlyEntities) {

		const int infoFrequency = _order.size() > 1000000 ? 100000 : 10000;
		vector<std::shared_ptr<CompactBCV>> bcvCache;
		bcvCache.resize(this->_weightedGraph->nodes.size(), 0);

		int counter = 0;
		for (vector<unsigned int>::iterator iter = _order.begin();
				iter < _order.end(); iter++) {
			const unsigned int focusWordGraphID = *iter;
			//		//only take specific one for debugging
			//		if (candidateNode.second != "<http://dbpedia.org/ontology/Province>"){
			//			continue;
			//		}

			shared_ptr<CompactBCV> combinedbcv = computeBCAIncludingEdgesCached(
					_weightedGraph, focusWordGraphID, bca_alpha, bca_eps,
					_predGraphIDs, bcvCache);

			//Note: this skip is done after the node has been added to the cache to make sure the value can be reused.
			if (onlyEntities) {
				if (!isEntity(this->_weightedGraph, focusWordGraphID)) {
					continue;
				}
			}

			const int focusWordGloveID = graphIDToGloveID(focusWordGraphID);

			for (vector<pair<unsigned int, double>>::const_iterator iter =
					combinedbcv->values.cbegin();
					iter != combinedbcv->values.cend(); iter++) {
				int contextWordGraphID = iter->first;
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter->second;
				CREC crec = CREC(focusWordGloveID, contextWordGloveID, freq);
				writeCrec(&crec, glove_input_file_out);
			}

			counter++;
			if ((counter % infoFrequency) == 0) {
				cout << currentTime() << "Processed " << counter << "/"
						<< _order.size() << " BCV computations" << endl;
			}
		}

	}

	void writeVocabFileIncludingEdges(FILE *glove_vocab_file_out) {
		this->writeVocabFileNoEdges(glove_vocab_file_out);

		//still need to write all predicates to the vocab file
		for (vector<flyString>::iterator it = _predicateLabels.begin();
				it < _predicateLabels.end(); it++) {
			fprintf(glove_vocab_file_out, "%s nofr\n", it->get().c_str());
		}
	}

};

} //end namespace co_occurence_computer

namespace KGloVe {

void parametrizedRun(Parameters &param) {
	param.check();

	int outputFileNumber = 0;
	for (std::vector<std::tuple<std::string, bool, bool>>::const_iterator graphs =
			param.graphs.begin(); graphs != param.graphs.end(); graphs++) {
		string graphInputPath(std::get<0>(*graphs));

		string graphFilename_replaced_dots = boost::replace_all_copy(
				boost::replace_all_copy(graphInputPath, ".", "_"), "/", "_");

		string outputPrefix = param.outputPrefix;

		const bool removeLiterals = std::get<1>(*graphs);
		const bool saveVocabulary = std::get<2>(*graphs);

		co_occurence_computer::Co_occurenceComputer c(graphInputPath.c_str(),
				removeLiterals);
		if (saveVocabulary) {

			string vocab_out_path(outputPrefix);
			for (auto includeEdges : param.includeEdges) {
				string start = "glove_vocab";
				if (includeEdges) {
					vocab_out_path += (start + "_with_edges_"
							+ graphFilename_replaced_dots + ".txt");
				} else {
					vocab_out_path += (start + "_no_edges_"
							+ graphFilename_replaced_dots + ".txt");
				}

				FILE *glove_vocab_file_out = fopen(vocab_out_path.c_str(), "w");
				if (!glove_vocab_file_out) {
					string err = string("vocab file could not be opened ")
							+ vocab_out_path;
					cerr << err << endl;
					cerr << "does the output directory exist?" << endl;
					throw err;
				}
				if (includeEdges) {
					c.writeVocabFileIncludingEdges(glove_vocab_file_out);
				} else {
					c.writeVocabFileNoEdges(glove_vocab_file_out);
				}
				fclose(glove_vocab_file_out);

			}

		}

		string fileNamePART1 = graphFilename_replaced_dots + "-"
				+ (removeLiterals ? "no_literals" : "with_literals");

		for (std::vector<weigher::GraphWeigher*>::const_iterator weigher =
				param.weighers.cbegin(); weigher != param.weighers.cend();
				weigher++) {
			weigher::GraphWeigher *forwardWeigher = *weigher;
			c.reWeigh(*forwardWeigher);

			string fileNamePART2 = fileNamePART1 + "-forwardWeigher_"
					+ forwardWeigher->getName();
			for (std::vector<double>::const_iterator alphas =
					param.alphas.cbegin(); alphas != param.alphas.cend();
					alphas++) {
				double alpha = *alphas;
				for (std::vector<double>::const_iterator epss =
						param.epss.cbegin(); epss != param.epss.cend();
						epss++) {
					double eps = *epss;
					for (std::vector<bool>::const_iterator onlyEntitiess =
							param.onlyEntities.cbegin();
							onlyEntitiess != param.onlyEntities.cend();
							onlyEntitiess++) {
						bool onlyEntities = *onlyEntitiess;
						for (std::vector<bool>::const_iterator includeEdgess =
								param.includeEdges.cbegin();
								includeEdgess != param.includeEdges.end();
								includeEdgess++) {
							bool includeEdges = *includeEdgess;

							//creating filename
							string outputFileName =
									"glove_input_file-" + fileNamePART2
											+ "-alpha_"
											+ boost::lexical_cast<std::string>(
													alpha) + "-eps_"
											+ boost::lexical_cast<std::string>(
													eps)
											+ (onlyEntities ?
													"-onlyEntities_yes" :
													"-onlyEntities_no")
											+ (includeEdges ?
													"-edges_yes" : "-edges_no")
											+ ".bin";

							string outputPath = outputPrefix + "/"
									+ outputFileName;

							FILE *glove_input_file_out = fopen(
									outputPath.c_str(), "w");
							if (!glove_input_file_out) {
								throw "coocurence file could not be written to";
							}
							if (includeEdges) {
								c.computeFrequenciesIncludingEdges(alpha, eps,
										glove_input_file_out, onlyEntities);
							} else {
								c.computeFrequenciesNoEdges(alpha, eps,
										glove_input_file_out, onlyEntities);
							}
							fclose(glove_input_file_out);
							outputFileNumber++;
						}
					}
				}
			}
		}
	}
}

}			//end KGloVe
