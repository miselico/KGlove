
TODO
----
Main
* Add command line arguments parsing option

Graph:
* Add option to modify graph
	* Add option to add inverse relations to graph
	* Add option to add the reversed edged to graph (keeping edge label as is)
* Add option to save weighted graph


KGloVe:
* Make different modes accessible
	* Normal
	* With inverse relations
	* With reversed edges

BCA
* Can be further optimized by calling itself recursively when something is missing in the cache. Cycles have to be taken care of, though == check whether increasing locality (instead of followig chains) in BCA order improves performance
* To be evaluated: at each BCV computation a bit of page rank is dropped. When reused, each time this bit is multiplied. Does this mean that in the end a lot of pagerank is lost?

Random walks (biased RDF2vec)
* Add generation of random walks and saving to file

C++ style issues
* check constness of iterators wherever possible
* make unsigned when it should  be.
* Currently C style file handling is used. Move to C++ style




