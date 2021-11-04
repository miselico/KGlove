

This repository contains code related to embedding entities in RDF graphs into a vector space.




Overview of files:
=====================


Main
------
Processes coommand line arguments, and calls the other parts



BCA
----
Bookmark coloring algorithm, and pushed-Bookmark coloring algorithm.
Implementation of approximate personal pagerank, including the option to include a cache, which is used for all pairs PPR computation.

GraphWeigher
-------------
Converters adding weights to a directed, labeled graph.
Other conversions could be done by these 'weighers' as well. 

MyMaxPriorityQueue
--------------------
Inspired by the priorityQueueu in glib. This priority queue supports changing of the priority of any element.
This is used in BCA, specifically for determining the ordering of flowing paint.


nTriplesParser
---------------
A quick 'n dirty parser for n-triples. No error handling, very sensitive to errors in input.
**Make sure the input is proper n-triples before attempting to use!**


RDF2Co_occurence
-------------------
The code in this unit creates the co-occurence matrices from RDF graphs. 
This is not yet completely integrated with glove, meaning that glove needs to be called on the output of this matrix.


Utils
-----




RDF2Walk -- curently not included
-----------
Helper functions for going from a n-triples file to a file containing walks


GraphWalker -- not currently included
-------------
Strategies for walking over a graph where edges have a weight

RandomWalkExperiments -- not currently included
---------------------
Experiments the generation of random walks



