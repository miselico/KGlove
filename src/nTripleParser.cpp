/*
 * ntripleParser.cpp
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 *
 *
 *
 *
 */

#include <assert.h>
#include <utility>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <memory>

#include "nTripleParser.h"
#include "utils.h"

using namespace std;
using namespace boost;

namespace {

bool DEBUGMODE = false;

class Triple: public std::tuple<const string, const string, const string> {
public:
	Triple(const string& S, const string& P, const string& O) :
			std::tuple<const string, const string, const string>(S, P, O) {
	}
	const string S() {
		return std::get<0>(*this);
	}
	const string P() {
		return std::get<1>(*this);
	}
	const string O() {
		return std::get<2>(*this);
	}
};

// A blank node is at the start of this line, parses it of and returns the remaining of the line.
pair<string, string> parseBlankNodeOf(const string & line) {
	std::size_t firstUnderScoreColon = line.find("_:", 0);

	assert(firstUnderScoreColon != string::npos && "_: expected not found ");
	string atBNStart = line.substr(firstUnderScoreColon);
	//We assume it is ended by space.

	int cutpoint = atBNStart.find(' ');
	string blankNode = atBNStart.substr(0, cutpoint);
	string rest = atBNStart.substr(cutpoint + 1);
	return pair<string, string>(blankNode, rest);
}

pair<string, string> parseSubjectOf(const string & startline) {
	int firstUnderScoreColon = startline.find("_:", 0);
	int firstAngular = startline.find('<', 0);

	if (firstUnderScoreColon != -1 && firstUnderScoreColon < firstAngular) {
		return parseBlankNodeOf(startline);
	}
	//otherwise it is a resource
	string line = startline.substr(startline.find('<'));
	int cutpoint = line.find('>');
	string S = line.substr(0, cutpoint + 1);
	return pair<string, string>(S, line.substr(cutpoint + 1));
}

string parseObject(const string & startline) {

	//the line now either contains one more resource or a literal
	//Note, there could be a '<' in the string and according to production of N3 also a " in resource

	int a = startline.find_first_not_of(' ');
	string line = startline.substr(a);

	if (line[0] == '_') {
		return parseBlankNodeOf(line).first;
	} else if (line[0] == '<') {
		line = line.substr(0, line.rfind('>') + 1);
		return line;
	} else if (line[0] == '"') {
		//literal
		if (DEBUGMODE && ((line.find("^^") != string::npos) || (line.find("@") != string::npos))) {
			cerr << "TODO language tags and datatypes are not supported properly, assumed as part of the literal : " << line.c_str() << endl;
		}
		//cut of final dot and spaces
		while (line[line.size() - 1] == ' ' || line[line.size() - 1] == '.') {
			line = line.substr(0, line.size() - 2);
		}
		return line;
	} else {
		cerr << " Parsing error, invalid object " << line.c_str() << endl;
		exit(1);
	}

}

Triple parsetripleLine(const string & startline) {
	//line = <http://www.wikidata.org/ontology#gcPrecision> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#DatatypeProperty> .
	//line = _:node1ap4j3o3kx4 <http://example.com> <http://example.com#hasname> "Test"

	pair<string, string> S_rest = parseSubjectOf(startline);

	string S = S_rest.first;

	//predicate is easy
	string line = S_rest.second.substr(S_rest.second.find('<'));
	string P = line.substr(0, line.find('>') + 1);
	line = line.substr(line.find('>') + 1);

	string O = parseObject(line);

	return Triple(S, P, O);
}

pair<std::shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > buildRDFGraphInternal(const string & filename, bool removeLiteral) {

//The graph
	std::shared_ptr<QuickGraph::LabeledGraph> graph(new QuickGraph::LabeledGraph);

	vector<QuickGraph::Node> & nodes = graph->nodes;

	nodes.reserve(10000);

//temporary keep track of all the nodes
	std::unordered_map<string, unsigned int> addedNodes;

	ifstream infile(filename);

	if(!infile.is_open())
	{
	  // error! maybe the file doesn't exist.
		cerr << "Input file " << filename << " not found, exiting!!" << endl;
		exit(7);
	}

	string line;
	int count = 0;

	//the flyweights will take care of multiple edges having the smae label.
	while (std::getline(infile, line)) {
		if (trim_copy(line).empty()) { //skip whitespace
			continue;
		}
		if (line.find('#', 0) == 0) {
			//comment
			continue;
		}
		count++;
		if (count % 1000000 == 0) {
			cout << currentTime() << "Read " << count << " triples" << endl;
		}
		Triple values = parsetripleLine(line);
		if (removeLiteral && values.O()[0] == '"') {
			//do not add this literal
			continue;
		}

		//we prune edges defined twice, if either subject or object are new, then it cannot be a double edge.
		bool twiceDefinedEdgePossible = true;

		int subjectIndex;
		{ //scoping for name clashes
			string subject = values.S();
			auto resS = addedNodes.find(subject);
			if (resS == addedNodes.end()) {
				subjectIndex = nodes.size();
				nodes.push_back(subject);
				addedNodes[subject] = subjectIndex;
				twiceDefinedEdgePossible = false;
			} else {
				subjectIndex = resS->second;
			}
		}
		unsigned int objectIndex;
		{ //scoping for name clashes
			string object = values.O();
			auto resO = addedNodes.find(object);

			if (resO == addedNodes.end()) {
				objectIndex = nodes.size();
				nodes.emplace_back(object);
				addedNodes[object] = objectIndex;
				twiceDefinedEdgePossible = false;
			} else {
				objectIndex = resO->second;
			}
		}
		//add ede if not second time it is defined
		vector<QuickGraph::Edge> & edges = nodes[subjectIndex].edges;
		bool twiceDefinedEdge = false;
		if (twiceDefinedEdgePossible) {
			for (vector<QuickGraph::Edge>::iterator edge = edges.begin(); edge < edges.end(); edge++) {
				if (edge->targetIndex == objectIndex) {
					if (edge->label == values.P()) {
						twiceDefinedEdge = true;
						break;
					}
				}
			}
		}
		if (!twiceDefinedEdge) {
			edges.emplace_back(values.P(), 0.0, objectIndex);
		}

	}

	cout << currentTime() << "Read " << count << " triples Altogether" << endl;

	infile.close();
	graph->pack();
	return pair<std::shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> >(graph, addedNodes);
}

}

namespace n3parser {

pair<std::shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > buildRDFGraph(const string & filename) {
	return buildRDFGraphInternal(filename, false);
}

pair<std::shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > buildRDFGraph(const string & filename, const bool removeLiterals) {
	return buildRDFGraphInternal(filename, removeLiterals);
}



pair<std::shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > buildRDFGraphIgnoreLiterals(const string & filename) {
	return buildRDFGraphInternal(filename, true);

}
}
