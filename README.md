# Minimal_Congestion_Spanning_Tree

------------------------------------------- General Description ---------------------------------------------------------

We include three algorithms (written in python) in the present version that give estimations to the minimal STC problem.

First algorithm sCD_p (or sCD_p_q):

A descent algorithm that provides estimations in polynomial time to the minimal STC problems in undirected simple graphs with edge weights (unweighted graphs are identified with graphs with constant weight = 1)
Arguments: (G,T,p)
G = graph, it is instanced as a list of edges. An edge is instanced as a list [u,v,w] where u,v are the adjacent vertices and w is the weight
T = An arbitrary spanning tree of G instanced as a list of edges in T
p = positive integer or np.inf, a parameter for the descent functional (in this case the Lp congestion as it is defined in the preprint "Minimal Lp congestion spanning trees on weighted graphs" by A. Castejón Lafuente, E. Estévez, C. Meniño Cotón and M.C. Somoza)

Second and third algorithms (LOC_BFS and ROC):

These are algorithms that search a good spanning tree in the dual graph of a planar graph. Edge congestions in planar graphs are computed from lenghts of cycles in the dual graph, this allows a local control of the Lp congestions of a spanning tree and allows recursive methods to construct spanning trees in the dual graph that are locally optimal for the congestion problem, these should be good estimations to the STC problem in typical graphs.
Arguments: (G,V,p)
G = Planar graph, in the present version we assume that edges are in correspondence with straight segments of the plane, i.e., the planar graph is geodesic
V = list of vertices in R^2 given as np.arrays
p =  positive integer or np.inf, a parameter for the local optimization functional (in this case the Lp congestion as it is defined in the preprint "Minimal Lp congestion spanning trees on weighted graphs" by A. Castejón Lafuente, E. Estévez, C. Meniño Cotón and M.C. Somoza)

The output of every algorithm is: the Lp congestion of a spanning tree, a list with the edge congestions of every edge, the spanning tree (given as a list of edges)

These three algorithms are described in the work "Minimal Lp congestion spanning trees on weighted graphs" by A. Castejón Lafuente, E. Estévez, C. Meniño Cotón and M.C. Somoza) that can be found in 

Recall also that p=1 gives spanning trees that are good estimations for the LSST problem (low stretch spanning tree) for unweighted graphs.

-------------------------------------------- Instruccions from python console -------------------------------------------------------------------------

Set the variable number_cores in the file Congestion_Git.py as you wish. This sets the number of threads to be used by multiprocess.

Import in python console mandatory libraries and the repository files: spanning_trees.py, test_graphs.py, Congestion_Git.py

If you want to test the sCD_p algorithm for the classical congestion in the hypercube graph H_7 beginning from a random spanning tree:

>>> G = test_graphs.hypercube_graph(7)
>>> T = spanning_trees.random_spanning_tree(G)
>>> stc, edge_congestion_list, output_spanning_tree = Congestion_Git.sCD_p(G,T,np.inf)

And just wait. You should see the progress printed in the console, for H_7 the computation should end in 5 to 10 minutes (in a mid range pc). For H_10 it may take a full day. It is recommended to make the first trys in small graphs to test the power or your pc.

You can also use the experimental algorithm sCD_p_q, this algorithm uses the Lp congestion as descending functional but keeps the best Lq congestion along the explored trees. This is an alternative way to estimate the classical STC problem setting q = np.inf and p large (usually p = 10 is good).

For planar graphs, you can also use the LOC_BFS and ROC algorithms.

If you want to use the ROC algorithm in a rectangular grid with 30x15 vertices in order to estimate the L1 congestion:

>>> G, V = test_graphs.grid_graph(30,15)
>>> stc, edge_congestion_list, output_spanning_tree = Congestion_Git.ROC(G,V,1)

In this case, the progress gives what polygon is used as base vertex in the dual graph for the recursive method and the current value of the Lp congestion. LOC_BFS and ROC are faster than sCD_p but in general give worse estimates for the Lp congestions (but they are reasonably good and perform better in euclidean tilings).

--------------------------------------------------- Warnings ---------------------------------------------------------------

The complexity of the algorithm sCD_p for p=np.infty and unweighted graphs is roughly upper estimated by O(m^3n^4) where m is the number of edges and n is the number of vertices.

For typicial graphs computer expriments give an effective complexity much slower (between O(m) and O(m^2)). Moreover the descent algorithm gets, in general, good estimations in time O(mn). Almost all the work load belongs to the latest steps of the descent instances. The algorithm prints the current congestion of the tree at each descent instance, it is a way to check the progress of the algorithm.

The previous discussion implies that these algorithms, in their present state, are feasible for mid range applications (graphs with edges in the range [0,10^4]).

The sCD_p algorithm, which is the most important one, uses the library multiprocess in order to parallelize some computations. The number of threads must be initialized in the file Congestion.py. The deafault value is 16.

Multiprocess is konwn to have some issues in Windows systems. If this is your case, then you should virtualize or dual boot a linux distribution.

Mandatory python libraries: numpy, multiprocess, random, itertools, matplotlib, pylab

--------------------------------------------------- To be done in near future ------------------------------------------------------

Keep in mind that this project is handled entirely by myself (C. Meniño Cotón). Although these algorithms were tested in a large family of graphs, bugs could appear.

Remark: I am not claiming that this is the fastest implementation possible and, in fact, I detected several points where the algorithms can be greatly improved. This will be done in the future and depending in the community response (and, of course, the available time).

Remark 2: test_graphs.py include a reasonable amount of test graphs to work with: Complete, complete multipartite, random, hypercube, discrete tori, cubic grid, rectangular grid, triangular grid (several configurations), hexagonal grid (several configurations), random planar and some weighted versions of these graphs. These families were chosen since the classical STC problem was already studied in most of them. These families will be also increased in the near future. In any case, any graph can be instanced by the user at any moment as a list of edges. Remark also that, in this program, vertices are labelled by intergers from 0 to n-1.

Remark 3: It is also expected to implement a graphical front end using PyQt.





