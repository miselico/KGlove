/*
 * GraphWeigher.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWeigher.h"

#include <iostream>

#include <string>
#include <unordered_map>
#include <fstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "graph/LabeledGraph.h"
#include "utils.h"
using namespace std;

typedef boost::flyweight<std::string> flyString;

namespace { //helpers

//count freq of each property:
unordered_map<flyString, double> absolute_predicate_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
	unordered_map<flyString, double> absolute_freq;
	std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	for (auto nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
		for (auto edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			flyString predictae = edgeI->label;
			auto oldValI = absolute_freq.find(predictae);
			if (oldValI == absolute_freq.end()) {
				absolute_freq[predictae] = 1;
			} else {
				absolute_freq[predictae] = oldValI->second + 1;
			}
		}
	}
	return absolute_freq;
}

////counts the absolute frequence of a the objects
//unordered_map<flyString, double> absolute_object_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
//	unordered_map<flyString, double> absolute_freq;
//	std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
//	for (auto nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
//		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
//		for (auto edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
//			flyString object = edgeI->target->label;
//			auto oldValI = absolute_freq.find(object);
//			if (oldValI == absolute_freq.end()) {
//				absolute_freq[object] = 1;
//			} else {
//				absolute_freq[object] = oldValI->second + 1;
//			}
//		}
//	}
//	return absolute_freq;
//}

//counts the absolute frequence of a the objects - adapted form the above as it seems more reasonable to keep an iterator, but then iteratrs cannot be hashed -> index
vector<double> absolute_object_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
	vector<double> absolute_freq;
	absolute_freq.resize(baseGraph->nodes.size(), 0);
	std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	for (auto nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
		for (auto edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			int objectIndex = edgeI->targetIndex;
			absolute_freq[objectIndex]++;
		}
	}
	return absolute_freq;
}

unordered_map<pair<flyString, flyString>, double, pairhash> absolute_predicate_object_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
	unordered_map<pair<flyString, flyString>, double, pairhash> absolute_freq;
	std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	for (auto nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
		for (auto edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			flyString predicate = edgeI->label;
			flyString object = nodes[edgeI->targetIndex].label;
			pair<flyString, flyString> pair(predicate, object);
			auto oldValI = absolute_freq.find(pair);
			if (oldValI == absolute_freq.end()) {
				absolute_freq[pair] = 1;
			} else {
				absolute_freq[pair] = oldValI->second + 1;
			}
		}
	}
	return absolute_freq;
}

//inverse all freq
template<class type, class Hash = std::hash<type>> void inverse_the_frequency_map(unordered_map<type, double, Hash> & absolute_freq) {
	for (typename unordered_map<type, double, Hash>::iterator iter = absolute_freq.begin(); iter != absolute_freq.end(); iter++) {
		absolute_freq[iter->first] = 1.0 / iter->second;
	}
}

//inverse all freq
void inverse_the_frequency_vec(vector<double> & absolute_freq) {
	for (typename vector<double>::iterator iter = absolute_freq.begin(); iter != absolute_freq.end(); iter++) {
		*iter = 1.0 / *iter;
	}
}

//normalize all weights in unbalanced (sum weight on outedges == 1)
void normalize(std::shared_ptr<QuickGraph::LabeledGraph> unbalanced) {

	for (vector<QuickGraph::Node>::iterator nodeI = unbalanced->nodes.begin(); nodeI != unbalanced->nodes.end(); nodeI++) {
		//normalize outgoing from this node
		double totalWeight = 0;
		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
		for (vector<QuickGraph::Edge>::const_iterator edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			totalWeight += edgeI->weight;
		}
		double totalWeightInverse = 1.0 / totalWeight;
		for (vector<QuickGraph::Edge>::iterator edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			double normalizedWeight = edgeI->weight * totalWeightInverse;
			edgeI->weight = normalizedWeight;
		}
	}

}

void setValueToOne(QuickGraph::LabeledGraph & /* baseGraph */, QuickGraph::Edge & edge) {
	edge.weight = 1;
}

} //and anonymous namespace helpers




namespace weigher{


unordered_map<string, double> readDBPediaPageRanks(string tsvFile) {
	unordered_map<string, double> ranks;

	ifstream infile(tsvFile);

	if (!infile.is_open()) {
		// error! maybe the file doesn't exist.
		cerr << "Input file " << tsvFile << " not found, exiting!!" << endl;
		exit(7);
	}

	string line;

	while (std::getline(infile, line)) {
		if (line.find_first_not_of(' ') == std::string::npos) {
			continue;
		}
		if (line.find_first_of('#') == 0) {
			//comment
			continue;
		}

		vector<string> SplitVec; // #2: Search for tokens
		boost::split(SplitVec, line, boost::is_any_of("\t"), boost::token_compress_off);
		string resource = "<" + SplitVec[0] + ">";
		double rank = boost::lexical_cast<double>(SplitVec[1]);

		ranks[resource] = rank;
	}
	infile.close();
	return ranks;
}


void UniformWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {

	baseGraph->forEachEdge(setValueToOne);
	normalize(baseGraph);

}

namespace {
class PredicateToWeight {
	const unordered_map<flyString, double> weights;
public:

	PredicateToWeight(const unordered_map<flyString, double> & weights) :
			weights(weights) {

	}
	void operator()(QuickGraph::LabeledGraph & /*baseGraph*/, QuickGraph::Edge & edge) {
		double newWeight = this->weights.at(edge.label);
		edge.weight = newWeight;
	}

};
}

void PredicateFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {

	unordered_map<flyString, double> absolute_freq = absolute_predicate_freq(baseGraph);
	baseGraph->forEachEdge(PredicateToWeight(absolute_freq));
	normalize(baseGraph);

}

void InversePredicateFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {

//count freq of each property:
	unordered_map<flyString, double> absolute_freq = absolute_predicate_freq(baseGraph);
	inverse_the_frequency_map(absolute_freq);
	baseGraph->forEachEdge(PredicateToWeight(absolute_freq));
	normalize(baseGraph);

}

void ObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	vector<double> absolute_freq = absolute_object_freq(baseGraph);
	PushDownWeigher subWeigher(absolute_freq);
	return subWeigher.weigh(baseGraph);
}

void InverseObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	vector<double> absolute_freq = absolute_object_freq(baseGraph);
	inverse_the_frequency_vec(absolute_freq);
	PushDownWeigher subWeigher(absolute_freq);
	return subWeigher.weigh(baseGraph);
}

void InverseObjectFrequencyWeigherSplitDown::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	vector<double> absolute_freq = absolute_object_freq(baseGraph);
	inverse_the_frequency_vec(absolute_freq);
	SplitDownWeigher subWeigher(absolute_freq);
	return subWeigher.weigh(baseGraph);
}

namespace {
void predObWeigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph, std::unordered_map<std::pair<flyString, flyString>, double, pairhash> & freqs) {
	std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	for (std::vector<QuickGraph::Node>::iterator nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
		for (std::vector<QuickGraph::Edge>::iterator edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			flyString pred = edgeI->label;
			flyString obj = nodes[edgeI->targetIndex].label;
			std::pair<flyString, flyString> pair(pred, obj);
			double weight = freqs.at(pair);
			edgeI->weight = weight;
		}
	}
	normalize(baseGraph);
}
}

void PredicateObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	std::unordered_map<std::pair<flyString, flyString>, double, pairhash> absolute_freq = absolute_predicate_object_freq(baseGraph);
	predObWeigh(baseGraph, absolute_freq);
}

void InversePredicateObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	std::unordered_map<std::pair<flyString, flyString>, double, pairhash> absolute_freq = absolute_predicate_object_freq(baseGraph);
	inverse_the_frequency_map(absolute_freq);
	predObWeigh(baseGraph, absolute_freq);
}

void PushDownWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	assert(nodes.size() == this->nodeWeights.size());

	for (unsigned int i = 0; i < this->nodeWeights.size(); i++) {
		std::vector<QuickGraph::Edge> & edges = nodes[i].edges;

		for (std::vector<QuickGraph::Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++) {

			int targetIndex = edge->targetIndex;
			double weight = this->nodeWeights[targetIndex];
			edge->weight = weight;
		}
	}
	normalize(baseGraph);
}

void PushDownWeigherMap::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	vector<QuickGraph::Node> & nodes = baseGraph->nodes;

	for (std::vector<QuickGraph::Node>::iterator iter = nodes.begin(); iter != nodes.end(); iter++) {
		QuickGraph::Node & theNode = *iter;
		std::unordered_map<string, double>::const_iterator findIter = this->nodeWeights.find(theNode.label);
		double weight;
		if (findIter == this->nodeWeights.cend()) {
			weight = this->defaultWeight;
		} else {
			weight = findIter->second;
		}
		std::vector<QuickGraph::Edge> & edges = theNode.edges;
		for (std::vector<QuickGraph::Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++) {
			edge->weight = weight;
		}
	}
	normalize(baseGraph);
}

void SplitDownWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	assert(nodes.size() == this->nodeWeights.size());

	//determine all in degrees
	vector<int> indegrees;
	indegrees.resize(nodes.size(), 0);
	for (unsigned int i = 0; i < this->nodeWeights.size(); i++) {
		std::vector<QuickGraph::Edge> & edges = nodes[i].edges;
		for (std::vector<QuickGraph::Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++) {
			int targetIndex = edge->targetIndex;
			indegrees[targetIndex]++;
		}
	}

	for (unsigned int i = 0; i < this->nodeWeights.size(); i++) {
		std::vector<QuickGraph::Edge> & edges = nodes[i].edges;

		for (std::vector<QuickGraph::Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++) {
			int targetIndex = edge->targetIndex;
			double weight = this->nodeWeights[targetIndex] / indegrees[targetIndex];
			edge->weight = weight;
		}
	}
	normalize(baseGraph);

}

}
