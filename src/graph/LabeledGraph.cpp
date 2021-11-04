/*
 * LabeledGraph.cpp
 *
 *  Created on: May 8, 2019
 *      Author: cochez
 */


#include "LabeledGraph.h"

namespace QuickGraph {

std::shared_ptr<QuickGraph::LabeledGraph> reverseGraph(std::shared_ptr<const QuickGraph::LabeledGraph> baseGraph) {
	std::shared_ptr<QuickGraph::LabeledGraph> newGraph(new QuickGraph::LabeledGraph);
	const std::vector<QuickGraph::Node> & originalNodes = baseGraph->nodes;
	std::vector<QuickGraph::Node> & newNodes = newGraph->nodes;

	//add all nodes:
	newNodes.reserve(originalNodes.size());
	for (std::vector<QuickGraph::Node>::const_iterator iter = originalNodes.begin(); iter != originalNodes.end(); iter++) {
		newNodes.emplace_back(iter->label);
	}

	//add all edges reversed

	for (std::vector<QuickGraph::Node>::const_iterator oldSrcNode = originalNodes.begin(); oldSrcNode != originalNodes.end(); oldSrcNode++) {
		const int oldSourceIndex = oldSrcNode - originalNodes.begin();
		const int newTargetIndex = oldSourceIndex;
		//const std::vector<QuickGraph::Node>::iterator newTargetNode = newNodes.begin() + newTargetIndex;

		std::vector<QuickGraph::Edge> oldEdges = oldSrcNode->edges;
		for (std::vector<QuickGraph::Edge>::const_iterator oldEdge = oldEdges.begin(); oldEdge != oldEdges.end(); oldEdge++) {
			int oldTargetIndex = oldEdge->targetIndex;
			int newSourceIndex = oldTargetIndex;
			QuickGraph::Node & newSourceNode = newNodes[newSourceIndex];
			newSourceNode.edges.emplace_back(oldEdge->label, oldEdge->weight, newTargetIndex);
		}
	}
	return newGraph;
}

}
