LoopFinder (Matheus Viana, vianamp@gmail.com - 14.05.15)
==========

Brute force loop detector in graphs. The user must specify the size L of the desired loop. This algorithm generates sequences of L interger numbers and test if these sequences correspond to closed loops in a given graph.

WARNING: Multiple connections, auto-loops and edges directionality are ignored in this algorithm.

INPUT
-----

Connection list in the format .gnet. In this format the first line specifies the total number of nodes and the subsequent lines correspond to the connections between nodes. The third column corresponds to the edge weight, which must be present but it is ignored in the current version.

[Example - sample.gnet]

12
0 2 1
1 3 1
2 3 1
2 4 1
3 4 1
4 5 1
5 6 1
5 9 1
6 7 1
6 8 1
6 9 1
8 9 1
8 11 1
9 10 1

HOW TO USE:
-----------

./LoopFinder filename_without_extension loop_size

[Example]

./LoopFinder sample 3

OUTPUT:
-------

The resulting analysis is shown on the screen and a xgmml file is saved in the same folder. On the screen it is shown the total number of loops found, what are the loop and the participation vector, which corresponds to the number of loops that each each node belongs to.

The xgmml file corresponds to the input network where the nodes attributes "width" and "height" are proportional to their participation. This file can be visualized using the software Cytoscape (http://www.cytoscape.org/), for instance.

[Example]

@Number of loops of length 3 found: 3.
@List of loops:
4 2 3 
9 5 6 
9 6 8 
@Participation of each node:
0
0
1
1
1
1
2
0
1
2
0
0