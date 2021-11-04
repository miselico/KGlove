/*
 * Main.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

//#include "RandomWalkExperiments.h"
#include "nTripleParser.h"

#include <iostream>

#include <utility>
#include <unordered_map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <fstream>
#include <boost/lexical_cast.hpp>

#include "RDF2Co_occurenceComplex.h"

using namespace std;

//int main(int argc, char **argv) {
//	//string filename = "wikidata-simple-statements-10_000000-sample.nt";
////	string filename = "../../datasets/dbPedia/allData27_30M.nt";
//	std::pair<std::shared_ptr<QuickGraph::LabeledGraph>, std::unordered_map<string, int>> g = n3parser::buildRDFGraph(filename);
//	cout << g.first->nodes.size() << endl;
//}

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

//int main(int argc, char **argv) {
//	try {
//		//char const * fileName = "dbpedia_base64_mtr100_mte100-train.nt";
//		char const * fileName = "freebase_mtr100_mte100-test.nt";
//		if (argc > 1) {
//			fileName = argv[1];
//		}
//		cout << "running2" << endl;
//		pair<shared_ptr<QuickGraph::LabeledGraph>, unordered_map<string, unsigned int> > graphAndNodeIndex = n3parser::buildRDFGraph(fileName, true);
//		std::shared_ptr<QuickGraph::LabeledGraph> weightedGraph = graphAndNodeIndex.first;
//
//		UniformWeigher w;
//		w.weigh(weightedGraph);
//
//		std::vector<QuickGraph::Node> & nodes = weightedGraph->nodes;
//		for (std::vector<QuickGraph::Node>::iterator nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
//			boost::flyweight<std::string> subject = nodeI->label;
//			std::vector<QuickGraph::Edge> & edges = nodeI->edges;
//			for (std::vector<QuickGraph::Edge>::iterator edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
//				boost::flyweight<std::string> predicate = edgeI->label;
//				boost::flyweight<std::string> object = nodes[edgeI->targetIndex].label;
//				long double weight = edgeI->weight;
//				std::cout << subject << ' ' << predicate << ' ' << object << ' ' << weight << std::endl;
//			}
//		}
//
//	} catch (char const* str) {
//		cout << str << endl;
//		throw str;
//	}
//	return 0;
//}

int COMPLEXmain(int argc, char **argv) {
	try {
		//char const * fileName = "dbpedia_base64_mtr100_mte100-train.nt";
//		char const * fileName = "freebase_mtr100_mte100-test.nt";
		char const * fileName = "SmallTest4.nt";
		if (argc > 1) {
			fileName = argv[1];
		}

		RDF2CO::ParameterizedRun::Parameters p;
		//p.graphs.push_back(std::tuple<string, bool, bool>("368303ALL_MergedMultiline_no-empty-lines_sort-uniq_error-boxer.nt", false, true));
		p.graphs.push_back(std::tuple<string, bool, bool>(fileName, true, true));
		weigher::UniformWeigher w;
		//PushDownWeigherMap w(readDBPediaPageRanks("pagerank_en_2016-04.tsv"), 0.2);
		p.weighers.push_back(std::pair<weigher::GraphWeigher&, weigher::GraphWeigher&>(w, w));
		p.alphas.push_back(0.3);
		p.epss.push_back(0.00001);
		p.normalize.push_back(true);
		p.onlyEntities.push_back(false);
		RDF2CO::ParameterizedRun::parametrizedUltimateRun(p);
	} catch (char const* str) {
		cout << str << endl;
		throw str;
	}
	return 0;
}

//int mainRandomWalk(int argc, char **argv) {
//
//
//
//
//
//	//RandomWalkExperiments::performExperiments(1);
//
//	//return 0;
//
//
//
//
//
//
//	if (argc < 2) {
//		cerr << "Usage: commandName [number], where number from" << endl;
//		cerr << "	1. Uniform" << endl << "	2. Predicate freq" << endl << "	3. Predicate inv. freq" << endl << "	4. Object freq" << endl << "	5. Object inv. freq" << endl << "	6. P-O freq" << endl
//				<< "	7. P-O inv. freq" << endl << "	8. Pagerank weight" << endl << "	9. inv Pagerank weight" << endl << "	10. Pagerank split weight" << endl << "	11. inv Pagerank split weight" << endl
//				<< "	12. Object inv. freq split" << endl << "	Object freq split == uniform"
//				<< endl;
//		return 1;
//	}
//
//	TStr strategy(argv[1]);
//	int strategyNumber = strategy.GetInt();
//
//	char * outFileName = NULL;
//	if (argc > 2){
//		outFileName = argv[2];
//	}
//
//	RandomWalkExperiments::performExperiments(strategyNumber, outFileName);
//
//
//
//	//RDF2CO::performExperiments();
//	return 0;
//}

