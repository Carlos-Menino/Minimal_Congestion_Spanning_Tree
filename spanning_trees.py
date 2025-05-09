#Basic methods to operate with spanning trees
#A graph is instanced as a list of edges. An edge is instanced as a list [u, v, w] where u is the initial vertex, v is the ending vertex and w is a positive weight
#Unweighted graphs are treated as graphs with constant weight w = 1
#At this moment we are not dealing with oriented graphs, thus we are not going to make distinction between [u,v,w] and [v,u,w].
# Please, if you use some part of this program, cite as:
# Meniño Cotón, C. (2025). Minimal Spanning tree estimation via descent methods (Version 1.0) [Computer software]. https://github.com/Carlos-Menino/Minimal_Congestion_Spanning_Tree
# Thanks!
import random as rd

def vertex_degree(G,v): #deg(v) for v vertex in G
    value = 0
    for edge in G:
        if edge[0] == v or edge[1] == v:
            value = value + 1
    return value

def weighted_degree(G,v): # sum of weights of adjacent edges to a vertex v
    value = 0
    for edge in G:
        if edge[0] == v or edge[1] == v:
            value = value + edge[2]
    return value

def vertex_star(G,v): #List of adjacent vertices to v
    star = []
    for edge in G:
        if edge[0] == v or edge[1] == v:
            star.append(edge)
    return star   

def vertex_star_vertex(G,v): #List of vertices connected to v by a single edge
    star = []
    for edge in G:
        if edge[0] == v:
            star.append(edge[1])
        elif edge[1] == v:
            star.append(edge[0])
    return star 

def graph_as_dictionary(G): #This method allows to change a graph to a dictionary format: G[v] is a list with triples [a,b,c], a is adjacent vertex to v, b is the edge-weight and c=+-1 depending in the orientation of the edge
    n = cardinal(G)
    Dic_G = dict()
    for i in range(n):
        Dic_G[i] = []
    for e in G:
        Dic_G[e[0]].append([e[1],e[2],1])
        Dic_G[e[1]].append([e[0],e[2],-1])
    return Dic_G

def cardinal(G): # Gets the number of vertices of a graph defined as a list
    V = []
    cardinal = 0
    for edge in G:
        if not(edge[0] in V):
            V.append(edge[0])
            cardinal = cardinal + 1
        if not(edge[1] in V):
            V.append(edge[1])
            cardinal = cardinal + 1
    return cardinal 

def vertex_set(Graph): # Gets the vertex set of a graph defined as a list
    V = []
    for e in Graph:
        if not(e[0] in V):
            V.append(e[0])
        if not(e[1] in V):
            V.append(e[1])
    return V

def labeling_tree(T): #assigns a tree orientation to a tree T
    T_copy = T.copy()
    T_label = [T[0]]
    V_label = [T[0][0],T[0][1]]
    V_level = [(T[0][0],0),(T[0][1],1)]
    T_copy.remove(T[0])
    while len(T_copy)!=0:
        for e in T_copy:
            if e[0] in V_label:
                T_label.append(e)
                V_label.append(e[1])
                index = V_label.index(e[0])
                V_level.append((e[1],V_level[index][1]+1))
                T_copy.remove(e)
                break
            elif e[1] in V_label:
                T_label.append([e[1],e[0],e[2]])
                V_label.append(e[0])
                index = V_label.index(e[1])
                V_level.append((e[0],V_level[index][1]+1))
                T_copy.remove(e)
                break
    #root = T[0][0]
    return T_label, V_label, V_level

def spanning_tree(G): # Generates a spanning tree T for a graph G
    T = [G[0]]
    VT = [T[0][0],T[0][1]]
    vertices = cardinal(G)
    G_copy = G.copy()
    G_copy.remove(G[0])
    while(len(VT)<vertices):
        edge = G_copy[0]
        if not(edge[0] in VT) and edge[1] in VT:
            T.append(edge)
            VT.append(edge[0])
            G_copy.remove(edge)
        if not(edge[1] in VT) and edge[0] in VT:
            T.append(edge)
            VT.append(edge[1])
            G_copy.remove(edge)
        else:
            G_aux = G_copy
            for i in range(len(G_aux)):
                G_copy[i-1] = G_aux[i]       
    return T

def random_spanning_tree(Graph): # Generates a random spanning tree T for a graph G
    Graph = rd.sample(Graph,len(Graph))    
    Tree = [Graph[0]]
    VT = [Tree[0][0],Tree[0][1]]
    vertices = cardinal(Graph)
    G_copy = Graph.copy()
    G_copy.remove(Graph[0])
    while(len(VT)<vertices):
        edge = G_copy[0]
        if not(edge[0] in VT) and edge[1] in VT:
            Tree.append(edge)
            VT.append(edge[0])
            G_copy.remove(edge)
        if not(edge[1] in VT) and edge[0] in VT:
            Tree.append(edge)
            VT.append(edge[1])
            G_copy.remove(edge)
        else:
            G_aux = G_copy.copy()
            for i in range(len(G_aux)):
                G_copy[i] = G_aux[i-1] 
    return Tree

def hypercube_spanning_tree(n): #Generates thee symmetric spanning tree of the hypercube
    T = [[0,1,1]]
    for i in range(1,n):
        T_copy = T.copy()
        for e in T_copy:
            T.append([e[0]+2**i,e[1]+2**i,1])
        T.append([0,2**i,1])
    return T

