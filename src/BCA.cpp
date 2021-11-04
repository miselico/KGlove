/*
 * BCA.cpp
 *
 *  Created on: Nov 24, 2016
 *      Author: cochez
 */

#include "BCA.h"
#include <unordered_set>

#include <iostream>
#include "MyMaxPriorityQueue.h"

using namespace std;

BCV::BCV(CompactBCV & from) {
	for (std::vector<std::pair<unsigned int, double>>::const_iterator iter = from.values.begin(); iter != from.values.begin(); iter++) {
		this->fixPaint(iter->first, iter->second);
	}
}

void BCV::fixPaint(unsigned int ID, double amount) {
	double startamount = 0;
	std::unordered_map<unsigned int, double>::iterator it = this->find(ID);
	if (it != this->end()) {
		startamount = it->second;
		double newAmount = startamount + amount;
		this->at(ID) = newAmount;
	} else {
#ifdef NDEBUG
		this->emplace(ID, amount);
#else
		auto ret = this->emplace(ID, amount);
		assert(ret.second);
#endif
	}
}

string BCV::toStringWithOnlyNodeLabels(const std::shared_ptr<QuickGraph::LabeledGraph> network) {
	string s = "{";
	string separator = "";

	for (unordered_map<unsigned int, double>::iterator iter = this->begin(); iter != this->end(); iter++) {
		s += separator;
		const int k = iter->first;
		string entity = network->nodes[k].label;
		s += entity;
		s += " = ";
		const double v = iter->second;
		s += v;
		separator = ", ";
	}

	s.append("}");
	return s;
}

void BCV::removeEntry(unsigned int ID) {
	this->erase(ID);
}

void BCV::normalizeInPlace() {
	double totalSum = 0.0;

	for (unordered_map<unsigned int, double>::iterator iter = this->begin(); iter != this->end(); iter++) {
		totalSum += iter->second;
	}
	for (unordered_map<unsigned int, double>::iterator iter = this->begin(); iter != this->end(); iter++) {
		double value = iter->second;
		double scaled = value / totalSum;
		this->emplace(iter->first, scaled);
	}
}

void BCV::add(const BCV & other) {
	for (unordered_map<unsigned int, double>::const_iterator iter = other.cbegin(); iter != other.cend(); iter++) {
		unsigned int ID = iter->first;
		double addition = iter->second;
		BCV::iterator ownValue = this->find(ID);
		if (ownValue == this->end()) {
			this->emplace(ID, addition);
		} else {
			this->emplace(ID, ownValue->second + addition);
		}

	}
}

void BCV::add(const CompactBCV & other) {
	for (vector<pair<unsigned int, double>>::const_iterator iter = other.values.cbegin(); iter != other.values.cend(); iter++) {
		unsigned int ID = iter->first;
		double addition = iter->second;
		BCV::iterator ownValue = this->find(ID);
		if (ownValue == this->end()) {
			this->emplace(ID, addition);
		} else {
			this->emplace(ID, ownValue->second + addition);
		}
	}
}

CompactBCV::CompactBCV(const BCV & bcv) {
	for (std::unordered_map<unsigned int, double>::const_iterator iter = bcv.cbegin(); iter != bcv.cend(); iter++) {
		values.emplace_back(iter->first, iter->second);
	}
	this->values.shrink_to_fit();
}


class BCAQueue: public MyMaxPriorityQueue<unsigned int> {
public:
	void addPaintTo(unsigned int toID, double paint) {
		double current = this->GetPriority(toID);
		this->SetPriority(toID, current + paint);
	}

	void addPaintToIfMoreAsEps(unsigned int toID, double paint, double eps) {
		double current = this->GetPriority(toID);
		double sum = current + paint;
		if (sum > eps) {
			this->SetPriority(toID, sum);
		}
	}

	bool empty() {
		return this->IsEmpty();
	}

	pair<unsigned int, double> pop() {
		double paint = this->GetMaxPriority();
		unsigned int ID = this->PopMax();
		return pair<unsigned int, double>(ID, paint);
	}

};

/**
 * Compute the bookmarking coloring algorithm (≃ personalized page rank) between node b_ID and all other nodes in the graph. Using teleportation parameter alpha and cut-off value eps.
 */
BCV computeBCA(std::shared_ptr<QuickGraph::LabeledGraph> network, unsigned int b_ID, double alpha, double eps) {
	BCAQueue Q;
	BCV p;
	Q.addPaintTo(b_ID, 1.0);
	while (!Q.empty()) {
		pair<int, double> element = Q.pop();
		int i = element.first;
		double w = element.second;
		p.fixPaint(i, alpha * w);
		if (w < eps) {
			continue;
		}
		std::vector<QuickGraph::Edge> & edges = network->nodes[i].edges;

		for (auto edge = edges.begin(); edge != edges.end(); edge++) {
			int j = edge->targetIndex;
			double edgeWeight = edge->weight;
			double paintToJ = (1.0 - alpha) * w * edgeWeight;
			Q.addPaintTo(j, paintToJ);
		}
	}
	return p;

}

/**
 * Compute the bookmarking coloring algorithm (≃ personalized page rank) between node b_ID and all other nodes in the graph. Using teleportation parameter alpha and cut-off value eps.
 *
 * This version also adds pagerank for the edges, each time when BCA traverses them.
 *
 * The first element of the pair is the normal BCV, the second one the pagerank for the predicates
 *
 */
BCV computeBCAIncludingEdges(std::shared_ptr<QuickGraph::LabeledGraph> network, unsigned int b_ID, double alpha, double eps, const std::unordered_map<flyString, unsigned int> & predIDs) {
	BCAQueue Q;
	BCV p;
	Q.addPaintTo(b_ID, 1.0);

	while (!Q.empty()) {
		pair<int, double> element = Q.pop();
		int i = element.first;
		double w = element.second;
		p.fixPaint(i, alpha * w);
		if (w < eps) {
			continue;
		}

		std::vector<QuickGraph::Edge> & edges = network->nodes[i].edges;

		for (auto edge = edges.begin(); edge != edges.end(); edge++) {
			int j = edge->targetIndex;
			double edgeWeight = edge->weight;
			double paintToJ = (1.0 - alpha) * w * edgeWeight;

			boost::flyweight<string> outEdgeLabel = edge->label;
			int outEdgeBCVID = predIDs.at(outEdgeLabel);
			//note: this was fixed. Not all paint going trough should be added to the edge. Only the fraction which will stay in the destination.
			p.fixPaint(outEdgeBCVID, paintToJ * alpha);
			//p.fixPaint(outEdgeBCVID, paintToJ);

			Q.addPaintTo(j, paintToJ);
		}
	}
	return p;
}

BCV computeBCACached(std::shared_ptr<QuickGraph::LabeledGraph> network, unsigned int b_ID, double alpha, double eps, std::unordered_map<unsigned int, BCV> & bcvCache) {
	BCAQueue Q;
	BCV p;
	Q.addPaintTo(b_ID, 1.0);
	while (!Q.empty()) {
		pair<int, double> element = Q.pop();
		int i = element.first;
		double w = element.second;
		std::unordered_map<unsigned int, BCV>::iterator precomputedI = bcvCache.find(i);
		if (precomputedI != bcvCache.end()) {
			BCV & precomputed = precomputedI->second;
			//TODO double check this:
			for (unordered_map<unsigned int, double>::iterator iter = precomputed.begin(); iter != precomputed.end(); iter++) {
				double scaled = iter->second * w;
				//Here the algorithm might have a slight difference with the original version.
				//The problem is that we cannot know the threshold directly because of the weights in the graph, this seems to be an okay estimation
				if (scaled > (eps * alpha)) {
					p.fixPaint(iter->first, scaled);
				}
			}
		} else {
			p.fixPaint(i, alpha * w);
			if (w < eps) {
				continue;
			}
			std::vector<QuickGraph::Edge> & edges = network->nodes[i].edges;

			for (auto edge = edges.begin(); edge != edges.end(); edge++) {
				int j = edge->targetIndex;
				double edgeWeight = edge->weight;
				double paintToJ = (1.0 - alpha) * w * edgeWeight;
				Q.addPaintTo(j, paintToJ);
			}
		}

	}
	bcvCache[b_ID] = p;
	return p;
}


void printBCV(BCV& bcv) {
	cout << "[";
	string sep = "";
	for (std::unordered_map<unsigned int, double>::const_iterator iter = bcv.cbegin(); iter != bcv.cend(); iter++) {
		cout << sep << iter->first << "," << iter->second;
		sep = " | ";
	}
	cout << "]" << endl;
}

shared_ptr<CompactBCV> computeBCAIncludingEdgesCached(std::shared_ptr<QuickGraph::LabeledGraph> network, unsigned int b_ID, double alpha, double eps,
		const std::unordered_map<flyString, unsigned int> & predIDs, vector<std::shared_ptr<CompactBCV>> & bcvCache) {
	BCAQueue Q;
	BCV p;
	Q.addPaintTo(b_ID, 1.0);

	while (!Q.empty()) {
		//printBCV(p);
		pair<unsigned int, double> element = Q.pop();
		unsigned int i = element.first;
		double w = element.second;

		std::shared_ptr<CompactBCV> cached = bcvCache[i];

		//TODO check whether non-filled elements indeed test false.
		if (cached) {
			CompactBCV * precomputed = cached.get();
			//TODO double check this:
			for (vector<pair<unsigned int, double>>::iterator iter = precomputed->values.begin(); iter != precomputed->values.end(); iter++) {
				double scaled = iter->second * w;
				//Here the algorithm might have a slight difference with the original version.
				//The problem is that we cannot know the threshold directly because of the weights in the graph, this seems to be an okay estimation
				if (scaled > (eps * alpha)) {
					p.fixPaint(iter->first, scaled);
				}
			}
		} else {
			p.fixPaint(i, alpha * w);
			if (w < eps) {
				continue;
			}
			std::vector<QuickGraph::Edge> & edges = network->nodes[i].edges;

			for (auto edge = edges.begin(); edge != edges.end(); edge++) {
				int j = edge->targetIndex;
				double edgeWeight = edge->weight;
				double paintToJ = (1.0 - alpha) * w * edgeWeight;

				flyString & outEdgeLabel = edge->label;
				int outEdgeBCVID = predIDs.at(outEdgeLabel);
				//note: this was fixed. Not all paint going trough should be added to the edge. Only the fraction which will stay in the destination.
				p.fixPaint(outEdgeBCVID, paintToJ * alpha);
				//p.fixPaint(outEdgeBCVID, paintToJ);
				Q.addPaintTo(j, paintToJ);

			}
		}
	}
	bcvCache[b_ID] = shared_ptr<CompactBCV>(new CompactBCV(p));				//   [b_ID] = p;
	return bcvCache[b_ID];
}

//From here on implementation of pushed Bookmark Coloring Algorithm

//void PBCV::fixPaint(pair<int, int> pred_obj_pair, double amount) {
//
//	double startamount = 0;
//	std::unordered_map<int, double>::iterator it = this->find(pred_obj_pair);
//	if (it != this->end()) {
//		startamount = it->second;
//	}
//
//	double newAmount = startamount + amount;
//	this->emplace(pred_obj_pair, newAmount);
//}
//
//PBCV computePBCA(std::shared_ptr<QuickGraph::LabeledGraph> network, int b_ID, double alpha, double eps) {
//
//	throw "not implemented";
//
////Note!! This queue contains the amount of paint which still has to be moved OUT of the nodes.
////Nothing of it should stay on the the node itself.
////This is the main confusing difference between the ordering of the algorithms.
//	PBCV p;
//	BCAQueue Q;
//	Q.addPaintTo(b_ID, 1.0);
//	while (!Q.empty()) {
//		pair<int, double> element = Q.pop();
//		int i = element.first;
//		double w = element.second;
//
////For all links i -> j
//		TNodeEdgeNet<string, WeightedPredicate>::TNodeI node_i = network->GetNI(i);
//		int node_i_outdeg = node_i.GetOutDeg();
//		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
//			int j = node_i.GetOutNId(outEdge);
//
//			WeightedPredicate edgeData = node_i.GetOutEDat(outEdge);
//			double edgeWeight = edgeData.W();
//			double paintToJ = w * edgeWeight;
//			cerr << "It seems the combinednode is not exactly what is should be! node_i.GetOutEId(outEdge) is not the edge label" << endl;
//			pair<int, int> combinedNode = pair<int, int>(node_i.GetOutEId(outEdge), j);
//			//stand-in for p_i = p_i + alpha*w
//			p.fixPaint(combinedNode, alpha * paintToJ);
//
//			double extrapaintToBeMovedOutFromJ = (1 - alpha) * paintToJ;
//			Q.addPaintToIfMoreAsEps(j, extrapaintToBeMovedOutFromJ, eps);
//		}
//	}
//	return p;
//}
