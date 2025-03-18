import random as rd

def vertex_degree(G,v):
    value = 0
    for edge in G:
        if edge[0] == v or edge[1] == v:
            value = value + 1
    return value

def weighted_degree(G,v):
    value = 0
    for edge in G:
        if edge[0] == v or edge[1] == v:
            value = value + edge[2]
    return value

def vertex_star(G,v):
    star = []
    for edge in G:
        if edge[0] == v or edge[1] == v:
            star.append(edge)
    return star   

def vertex_star_vertex(G,v):
    star = []
    for edge in G:
        if edge[0] == v:
            star.append(edge[1])
        elif edge[1] == v:
            star.append(edge[0])
    return star 

def graph_as_dictionary(G):
    n = cardinal(G)
    Dic_G = dict()
    for i in range(n):
        Dic_G[i] = []
    for e in G:
        Dic_G[e[0]].append([e[1],e[2],1])
        Dic_G[e[1]].append([e[0],e[2],-1])
    return Dic_G

def cardinal(G):
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

def vertex_set(Graph):
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

def global_tree(Span_Tree, SubTree):
    New_Tree = SubTree.copy()
    Vertices_SubTree = vertex_set(SubTree)
    Aux_Tree = Span_Tree.copy()
    while len(Vertices_SubTree)< len(Span_Tree)+1:
        e = Aux_Tree[0]
        if e[0] in Vertices_SubTree and not(e[1] in Vertices_SubTree):
            New_Tree.append(e)
            Vertices_SubTree.append(e[1])
            Aux_Tree.remove(e)
        if e[1] in Vertices_SubTree and not(e[0] in Vertices_SubTree):
            New_Tree.append(e)
            Vertices_SubTree.append(e[0])
            Aux_Tree.remove(e)
        else:
            Aux_copy = Aux_Tree.copy()
            for i in range(len(Aux_Tree)):
                Aux_Tree[i] = Aux_copy[i-1]
    return New_Tree

def spanning_tree(G): 
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

def spanning_subtree(T,A):
    TA = []
    for e in T:
        if (e[0] in A) and (e[1] in A):
            TA.append(e)
    return TA

def random_spanning_tree(Graph): #puede no generar un árbol tal y como está escrito [CORREGIR]
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

def hypercube_spanning_tree(n):
    T = [[0,1,1]]
    for i in range(1,n):
        T_copy = T.copy()
        for e in T_copy:
            T.append([e[0]+2**i,e[1]+2**i,1])
        T.append([0,2**i,1])
    return T

def good_spanning_tree(G): #A veces no sale árbol generador, revisar
    G = rd.sample(G,len(G))    
    V = []
    for edge in G:
        if not(edge[0] in V):
            V.append(edge[0])
        if not(edge[1] in V):
            V.append(edge[1])
    degrees = []
    for v in V:
        degrees.append(vertex_degree(G, v))
    d_max = max(degrees)
    max_index = degrees.index(d_max)
    root_vertex = V[max_index]
    VL=[root_vertex]
    init_set = [root_vertex]
    T = []
    while len(T)< len(V)-1:
        relative_degrees = []
        relative_stars = []
        for v in init_set:
            Star = vertex_star(G, v)
            aux = Star.copy()
            for e in Star:
                if e in T:
                    aux.remove(e)
                elif (e[0] == v) and (e[1] in VL):
                    aux.remove(e)
                elif (e[1] == v) and (e[0] in VL):
                    aux.remove(e)           
            relative_degrees.append(len(aux))
            relative_stars.append(aux)
        rel_max = max(relative_degrees)
        rel_index = relative_degrees.index(rel_max)
        for edge in relative_stars[rel_index]:
            T.append(edge)
            if edge[0] == init_set[rel_index]:
                VL.append(edge[1])
                init_set.append(edge[1])
            elif edge[1] == init_set[rel_index]:
                VL.append(edge[0])
                init_set.append(edge[0])
        init_set.remove(init_set[rel_index])
    return T

