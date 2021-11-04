/*
 * GraphWeigher.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef GRAPHWEIGHER_HA_
#define GRAPHWEIGHER_HA_

#include <memory>
#include "graph/LabeledGraph.h"
#include <unordered_map>

namespace weigher{

class GraphWeigher {

std::unordered_map<std::string, double> readDBPediaPageRanks(std::string tsvFile);

protected:
	GraphWeigher() {

	}

public:
	virtual ~GraphWeigher() {

	}
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const = 0;

	virtual std::string getName() const = 0;
};

class UniformWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;

	virtual std::string getName() const {
		return "UniformWeigher";
	}
};

class InversePredicateFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;

	virtual std::string getName() const {
		return "InversePredicateFrequencyWeigher";
	}
};

class PredicateFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
	virtual std::string getName() const {
		return "PredicateFrequencyWeigher";
	}
};

class ObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
	virtual std::string getName() const {
		return "ObjectFrequencyWeigher";
	}
};

class InverseObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
	virtual std::string getName() const {
		return "InverseObjectFrequencyWeigher";
	}
};

class PredicateObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
	virtual std::string getName() const {
		return "PredicateObjectFrequencyWeigher";
	}
};

class InversePredicateObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
	virtual std::string getName() const {
		return "InversePredicateObjectFrequencyWeigher";
	}
};

class InverseObjectFrequencyWeigherSplitDown: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
	virtual std::string getName() const {
		return "InverseObjectFrequencyWeigherSplitDown";
	}
};

/**
 * Assigns to each inedge the weight assigned to the nodes.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class PushDownWeigher: public GraphWeigher {
	const std::vector<double> nodeWeights;

public:
	/**
	 * All weights must be in the vector
	 */
	PushDownWeigher(const std::vector<double> & nodeWeights) :
			nodeWeights(nodeWeights) {
	}

	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;

	virtual std::string getName() const {
		return "Pushdownweigher";
	}
};

/**
 * Assigns to each inedge the weight assigned to the nodes.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class PushDownWeigherMap: public GraphWeigher {
	const std::unordered_map<std::string, double> nodeWeights;
	const double defaultWeight;

public:
	/**
	 * All weights must be in the vector
	 */
	PushDownWeigherMap(const std::unordered_map<std::string, double> & nodeWeights, double defaultWeight) :
			nodeWeights(nodeWeights), defaultWeight(defaultWeight) {
	}

	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;

	virtual std::string getName() const {
		return "PushdownweigherMap";
	}
};


/**
 * Assigns to each inedge the weight assigned to the nodes divided by the number of in edges.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight divided by #inedge
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node / #inedge
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class SplitDownWeigher: public GraphWeigher {
	const std::vector<double> nodeWeights;

public:
	/**
	 * All weights must be in the vector
	 */
	SplitDownWeigher(const std::vector<double> & nodeWeights) :
			nodeWeights(nodeWeights) {
	}

	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;

	virtual std::string getName() const {
		return "Splitdownweigher";
	}
};

}
#endif /* GRAPHWEIGHER_HA_ */
