/*
 * utils.h
 *
 *  Created on: Feb 18, 2017
 *      Author: cochez
 */

#ifndef UTILS_HA_
#define UTILS_HA_

#include <string>

static std::string RDF_TYPE("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static std::string OWL_THING("<http://www.w3.org/2002/07/owl#Thing>");

struct pairhash {
public:
	template<typename T, typename U>
	std::size_t operator()(const std::pair<T, U> &x) const {
		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};


//the following is taken from by http://stackoverflow.com/a/16358264

#include <iostream>
#include <ctime>

inline std::string currentTime() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, sizeof(buffer), "%Y-%m-%d %I:%M:%S ", timeinfo);
	std::string str(buffer);

	return str;

}




#endif /* UTILS_HA_ */
