# Please, if you use some part of this program, cite as:
# Meniño Cotón, C. (2025). Minimal Spanning tree estimation via descent methods (Version 1.0) [Computer software]. https://github.com/Carlos-Menino/Minimal_Congestion_Spanning_Tree
# Thanks!
import numpy as np
import random as rd
from multiprocess import Pool
import test_graphs
import spanning_trees
from itertools import repeat

number_cores = 16
        
def path_tree_vertices(T,V_label,V_level,v,w): #finds the unique path in a labeled spanning tree T joining v to w (v, w must be different)
     index_v = V_label.index(v)
     index_w = V_label.index(w)
     level_v = V_level[index_v][1]
     level_w = V_level[index_w][1]
     path = []
     ends = [T[i][1] for i in range(len(T))]
     signal = False
     v_aux = v
     w_aux = w
     level_v_aux = level_v
     level_w_aux = level_w
     while signal == False:
        if level_v_aux > level_w_aux:
             index_v_aux = ends.index(v_aux)
             path.append(T[index_v_aux])
             v_aux = T[index_v_aux][0]
             level_v_aux = level_v_aux - 1
        elif level_w_aux > level_v_aux:
             index_w_aux = ends.index(w_aux)
             path.append(T[index_w_aux])
             w_aux = T[index_w_aux][0]
             level_w_aux = level_w_aux - 1
        elif level_v_aux == level_w_aux:
            index_v_aux = ends.index(v_aux)
            path.append(T[index_v_aux])
            v_aux = T[index_v_aux][0]
            level_v_aux = level_v_aux - 1
            index_w_aux = ends.index(w_aux)
            path.append(T[index_w_aux])
            w_aux = T[index_w_aux][0]
            level_w_aux = level_w_aux - 1
        if v_aux == w_aux:
            signal = True
     return path          
                
def unweight(G):
    Gu = [[G[i][0],G[i][1]] for i in range(len(G))]
    return Gu

def subgraph(Graph, Vertex_Set):
    SubGraph = []
    for e in Graph:
        if (e[0] in Vertex_Set) and (e[1] in Vertex_Set):
            SubGraph.append(e)
    return SubGraph


def congestion_vector(G,T, VL, VN): # G weighted graph, T spanning tree, VL, VN --> tree labeling
    num_edges_tree = len(T)    
    congestion_vector = np.zeros(num_edges_tree).tolist()
    for e in G:
        if e in T:
            index = T.index(e)
            congestion_vector[index] += e[2]
        elif [e[1],e[0],e[2]] in T:
            index = T.index([e[1],e[0],e[2]])
            congestion_vector[index] += e[2]
        else:
         v = e[0]
         w = e[1]
         path = path_tree_vertices(T, VL, VN, v, w)
         for edge in path:
             index = T.index(edge)
             congestion_vector[index] += e[2]
    return congestion_vector

def L1_tree_congestion(G,T):
    New_T, New_VL, New_VN  = spanning_trees.labeling_tree(T)
    congestion_list = congestion_vector(G, New_T, New_VL, New_VN)
    L1_congestion = sum(congestion_list)
    return L1_congestion

def Linf_tree_congestion(G,T):
    New_T, New_VL, New_VN  = spanning_trees.labeling_tree(T)
    congestion_list = congestion_vector(G, New_T, New_VL, New_VN)
    Linf_congestion = max(congestion_list)
    return Linf_congestion

def process_path(G,T,tree_path,e,edge,p):
    new_edge_hurt_list = []  
    new_path_hurt_list = []
    index_edge = T.index(edge)
    New_T = T.copy()
    New_T[index_edge]=e #New tree obtained from T changing edge by e
    New_T, New_VL, New_VN  = spanning_trees.labeling_tree(New_T)
    new_edge_hurt_list = congestion_vector(G, New_T, New_VL, New_VN)
    
    aux_path = tree_path.copy()
    aux_path.remove(edge)
    aux_path.append(e)
    for aux_edge in aux_path:
        if aux_edge in New_T:
            aux_index = New_T.index(aux_edge)
        else:
            aux_index = New_T.index([aux_edge[1],aux_edge[0],aux_edge[2]])
        aux_hurt = new_edge_hurt_list[aux_index]
        new_path_hurt_list.append(aux_hurt)
        #new_edge_hurt_list[aux_index] = aux_hurt
    if p == np.inf :
        New_congestion = max(new_edge_hurt_list) 
        new_cycle_hurt = max(new_path_hurt_list)
    else:
        New_congestion = (sum([h**p for h in new_edge_hurt_list]))
        new_cycle_hurt = sum([h**p for h in new_path_hurt_list])
    return [New_congestion, new_cycle_hurt, new_edge_hurt_list, New_T, New_VL,New_VN]

#################### AUX FUNCTIONS ##############################################
def complete_multipartite_bound(order_list,p): #Gets an upperbound for the Lp, 1<= p < np.inf congestion in complete multipartite graphs K_{order_list}, conjecturally it is the Lp congestion
    suma = 0
    k = len(order_list)
    N = sum(order_list)
    if order_list[-1] == 1:
        for i in range(k-1):
            #for j in range(i+1,k-1):
                #suma = suma + 2*order_list[i]*order_list[j]
            suma = suma + order_list[i]*(N-order_list[i])**p
        #suma = suma + sum(order_list) - order_list[-1]
    else:
        for i in range(k):
            suma = suma + order_list[i]*(N-order_list[i])**p
        n1 = order_list[0]
        nk = order_list[-1]
        suma1 = suma + (nk-1)*((2*N-n1-nk-2)**p-(N-n1)**p) - (N-nk)**p  
        suma2 = suma - (N-n1)**p-(N-nk)**p + (nk*(N-nk-1)+2-n1)**p
        suma = min(suma1,suma2)
    return suma

#################### CONGESTION ALGORITHMS ###########################################

def sCD_p(G,T,p):
    T, VL,VN = spanning_trees.labeling_tree(T)
    print('Check OK')
    edge_hurt_list = congestion_vector(G, T, VL,VN)
    if p == np.inf:
        print('Computing classical congestion')
        current_congestion = max(edge_hurt_list) 
    else:
        print('Computing L',p, ' congestion')
        current_congestion = (sum([h**p for h in edge_hurt_list]))
    token = False
    while token == False:
        progress = 0
        count = 0
        check = False
        G = rd.sample(G,len(G))  
        for e in G:
            progress = progress + 1
            if not(e in T) and not([e[1],e[0],e[2]]in T):
                v = e[0]
                w = e[1]
                tree_path = path_tree_vertices(T, VL, VN, v, w)
                path_hurt_list = []
                for edge in tree_path: 
                    index_edge = T.index(edge)
                    path_hurt_list.append(edge_hurt_list[index_edge])
                    if p == np.inf:
                        current_cycle_hurt = max(path_hurt_list)
                    else:
                        current_cycle_hurt = sum([h**p for h in path_hurt_list])
                pool = Pool(number_cores)
                new_congestion_list = pool.starmap(process_path, zip(repeat(G),repeat(T),repeat(tree_path),repeat(e),tree_path,repeat(p)))
                path_congestions = [new_congestion_list[i][1] for i in range(len(new_congestion_list))]
                optimal_change = min(path_congestions)
                optimal_index = path_congestions.index(optimal_change)
                optimal_congestion = new_congestion_list[optimal_index][0]

                if (new_congestion_list[optimal_index][0] >= current_congestion):
                    count = count + 1
                else:    
                    current_congestion = optimal_congestion
                    New_T = new_congestion_list[optimal_index][3]
                    edge_hurt_list = new_congestion_list[optimal_index][2]
                    T, VL, VN = New_T, new_congestion_list[optimal_index][4], new_congestion_list[optimal_index][5]
                    print('Current L', p, ' congestion: ', current_congestion,'Current classical congestion:', max(edge_hurt_list))
                    count = 0
                    check = True
                        
                    break
                if (optimal_change < current_cycle_hurt)and(optimal_congestion >= current_congestion): #
                    current_congestion = optimal_congestion
                    New_T = new_congestion_list[optimal_index][3]
                    edge_hurt_list = new_congestion_list[optimal_index][2]
                    T, VL, VN = New_T, new_congestion_list[optimal_index][4], new_congestion_list[optimal_index][5]
                    print('Current L', p, ' congestion: ',current_congestion,  ', Current classical congestion:', max(edge_hurt_list)', edge count= ', count,)
                    count = 0
                    check = True
                   
                    break 
        if check == False:
            print('Current L', p, ' congestion: ',current_congestion, 'Current classical congestion:', max(edge_hurt_list))
            token = True
        if count == len(G)-len(T)+2:
            token = True
            print('Halt Stop')
    return current_congestion, edge_hurt_list, T


def sCD_p_q(G,T, p, q):
    T, VL,VN = spanning_trees.labeling_tree(T)
    print('Check OK')
    edge_hurt_list = congestion_vector(G, T, VL,VN)
    if q == np.inf:
        current_congestion = sum([h**p for h in edge_hurt_list])
        best_congestion = max(edge_hurt_list)
        best_T = T
        
    else:
        current_congestion = sum([h**p for h in edge_hurt_list])
        best_congestion = sum([h**q for h in edge_hurt_list])
        best_T = T
    token = False
    
    while token == False:
        progress = 0
        count = 0
        check = False
        G = rd.sample(G,len(G))  
        for e in G:
            progress = progress + 1
            if not(e in T) and not([e[1],e[0],e[2]]in T):
                v = e[0]
                w = e[1]
                tree_path = path_tree_vertices(T, VL, VN, v, w)
                path_hurt_list = []
                for edge in tree_path:
                    index_edge = T.index(edge)
                    path_hurt_list.append(edge_hurt_list[index_edge])
                    if p == 'inf':
                        current_cycle_hurt = max(path_hurt_list)
    
                    else:
                        current_cycle_hurt = sum([h**p for h in path_hurt_list])
                pool = Pool(number_cores)
                new_congestion_list = pool.starmap(process_path, zip(repeat(G),repeat(T),repeat(tree_path),repeat(e),tree_path,repeat(p)))
                path_congestions = [new_congestion_list[i][1] for i in range(len(new_congestion_list))]
                if q == 'inf':
                    congestion_list = [max(new_congestion_list[i][2]) for i in range(len(new_congestion_list))]
                    MIN = min(congestion_list)
                    if MIN < best_congestion:
                        best_congestion = MIN
                        best_T = new_congestion_list[congestion_list.index(MIN)][3]
                else:
                    congestion_list = [sum([h**q for h in new_congestion_list[i][2]]) for i in range(len(new_congestion_list))]
                    MIN = min(congestion_list)
                    if MIN < best_congestion:
                        best_congestion = MIN
                        best_T = new_congestion_list[congestion_list.index(MIN)][3]
                optimal_change = min(path_congestions)
                optimal_index = path_congestions.index(optimal_change)
    
                if (new_congestion_list[optimal_index][0] >= current_congestion):
                    count = count + 1
    
                else:
                    count = 0
                    current_congestion = new_congestion_list[optimal_index][0]
                    New_T = new_congestion_list[optimal_index][3]
                    new_edge_hurt_list = new_congestion_list[optimal_index][2]
                    T, VL, VN = spanning_trees.labeling_tree(New_T)
                    for arrow in New_T:
                        if arrow in T:
                            T_index = T.index(arrow)
                            New_T_index = New_T.index(arrow)
                            edge_hurt_list[T_index] = new_edge_hurt_list[New_T_index]
                        elif [arrow[1],arrow[0],arrow[2]]in T:
                            T_index = T.index([arrow[1],arrow[0],arrow[2]])
                            New_T_index = New_T.index(arrow)
                            edge_hurt_list[T_index] = new_edge_hurt_list[New_T_index]
                    if (q == np.inf): 
                        Linf = max(edge_hurt_list)
                        if best_congestion > Linf:
                            best_congestion = Linf
                            best_T = T
                        print('Best Linf congestion: ', best_congestion)
                    else:
                        Lq = sum([h**q for h in edge_hurt_list])
                        if best_congestion > Lq:
                            best_congestion = Lq
                            best_T = T
                        print('Best Lq congestion: ', best_congestion**(1/q))
                    check = True
                    
                    break 
                if (optimal_change < current_cycle_hurt):
                    current_congestion = new_congestion_list[optimal_index][0]
                    New_T = new_congestion_list[optimal_index][3]
                    new_edge_hurt_list = new_congestion_list[optimal_index][2]
                    T, VL, VN = spanning_trees.labeling_tree(New_T)
                    for arrow in New_T:
                        if arrow in T:
                            T_index = T.index(arrow)
                            New_T_index = New_T.index(arrow)
                            edge_hurt_list[T_index] = new_edge_hurt_list[New_T_index]
                        elif [arrow[1],arrow[0],arrow[2]]in T:
                            T_index = T.index([arrow[1],arrow[0],arrow[2]])
                            New_T_index = New_T.index(arrow)
                            edge_hurt_list[T_index] = new_edge_hurt_list[New_T_index]
                    if (q == np.inf): 
                        Linf = max(edge_hurt_list)
                        if best_congestion > Linf:
                            best_congestion = Linf
                            best_T = T
                        print('Best Linf congestion: ', best_congestion)
                    else:
                        Lq = sum([h**q for h in edge_hurt_list])
                        if best_congestion > Lq:
                            best_congestion = Lq
                            best_T = T
                        print('Best Lq congestion: ', best_congestion**(1/q))
                    check = True
                    count = 0
                    check = True
                   
                    break 
        if check == False:
            token = True
        if count == len(G)-len(T)+2:
            token = True
            print('Halt Stop')
    return best_congestion, edge_hurt_list, best_T

################## FOR PLANAR GRAPHS #############################################

def LOC_BFS(G,vertex_list,p): #Depurar programa
    DG, PL = test_graphs.dual_graph(G,vertex_list) #DG = dual graph, PL = cell list  
    last_congestion = np.inf #este lo usaremos para romper el while cuando aparezca una congestión peor
    ## Dual graph as dictionary
    Gd = spanning_trees.graph_as_dictionary(DG)
    LOC_congestion_lists = []
    for k in range(len(PL)):
        print('Polygon basis: ',k)
        Dual_Copy = DG.copy()
        S = [k]
        DT = []
        BS =[]
        for v in Gd[k]:
            if not(v[0]in S):
                S.append(v[0])
                BS.append(v[0])
                if v[2]==1:
                    DT.append([k,v[0],v[1]])
                else:
                    DT.append([v[0],k,v[1]])
            else:
                for e in DT: #Intentar que esta parte no sea un bucle, aunque en este punto no afecta a la velocidad
                    if (e[1]==v[0] or e[0] == v[0]) and e[2]>v[1]: #Check if this line works in a weighted graph
                        if v[2] == 1:
                            e = [k,v[0],v[1]]
                        else:
                            e = [v[0],k,v[1]]
        L = []
        for i in range(len(DG)):
            if (DG[i][1] in S) and ((DG[i][0] in S)) and not(DG[i] in DT): 
                L.append(DG[i])
            elif (DG[i] in DT)  and  i!=DG.index(DG[i]): #para contar aristas múltiples
                L.append(DG[i])     
        new_DT, VL,VN = spanning_trees.labeling_tree(DT)
        congestion_list = []
        for e in L:
            Dual_Copy.remove(e)
            v = e[0]
            w = e[1]
            tree_path = path_tree_vertices(new_DT, VL, VN, v, w)
            c = sum([edge[2] for edge in tree_path]) + e[2]
            congestion_list.append(c)
        Check_Out = False
        while ((len(DT)+1) < len(PL))and Check_Out == False:   ### AQUÍ EMPIEZA LO IMPORTANTE ###
            New_BS = []
            Join_Edges = dict()
            for v in BS:
                for w in Gd[v]:
                    if not(w[0] in S)and not(w[0]in New_BS):
                        New_BS.append(w[0])
                        Join_Edges[w[0]] = []
            for v in BS:
                for w in Gd[v]:
                    if w[0] in New_BS:
                        if w[2] == 1:
                            Join_Edges[w[0]].append([v,w[0],w[1]])
                        else:
                            Join_Edges[w[0]].append([w[0],v,w[1]])
            ## Possible BFS Extension
            B_Edges = []
            B_Edge_congestion = []
            
            New_BS = rd.sample(New_BS,len(New_BS)) # Hacemos una elección aleatoria del orden de vértices en la frontera para añadir aleatoriedad
            
            for w in New_BS:
                B_Edges.append(Join_Edges[w][0]) ##Primeras aristas de join edge realizan una BFS extensión del árbol
            ## Computing Extended Congestions
            New_DT = DT + B_Edges #BFS extensión
            
            BS = New_BS  # Actualizamos la frontera del árbol (no hace falta hacerlo aquí)
            S = S + New_BS # Actualizamos vértices en el árbol
            new_DT, VL,VN = spanning_trees.labeling_tree(New_DT)
            for w in New_BS:
                if  len(Join_Edges[w])>1: #switch point (si no es switch point no hay nada que hacer)
                    aux = []
                    for e in Join_Edges[w][1:len(Join_Edges[w])]:
                        v = e[0]
                        w = e[1]
                        tree_path = path_tree_vertices(new_DT, VL, VN, v, w)
                        c = sum([edge[2] for edge in tree_path]) + e[2]
                        aux.append(c)
                    B_Edge_congestion.append(aux) #Congestiones de aristas para la extensión anterior
            #Exterior Edges: (Conectan vértices de la frontera)
            L=[]
            ext_con = []
            for e in Dual_Copy:
                if (e[1] in New_BS) and ((e[0] in New_BS)): 
                    L.append(e)
            if len(L) != 0:
                for e in L:        
                    Dual_Copy.remove(e)
                    v = e[0]
                    w = e[1]
                    tree_path = path_tree_vertices(new_DT, VL, VN, v, w)
                    c = sum([edge[2] for edge in tree_path]) + e[2]
                    ext_con.append(c)
            
            if p < np.inf: #Cálculo de congestiones relativas del árbol extendido
                if len(B_Edge_congestion) == 0:
                    con = sum([c**p for c in (congestion_list + [0])]) + sum([c**p for c in (ext_con+[0])])
                    if con != 0:
                        rel_con = con/(len(congestion_list) + len(ext_con))
                else:    
                    con = sum([c**p for c in (congestion_list + [0])]) + sum([c**p for c in (ext_con+[0])]) + sum([sum([c**p for c in (E+[0])]) for E in B_Edge_congestion])
                    rel_con = con/(len(congestion_list) + sum([len(E) for E in B_Edge_congestion])+ len(ext_con))
            else:
                if len(B_Edge_congestion) == 0:
                    con = max(max(congestion_list + [0]), max(ext_con+[0]))
                    rel_con = con
                else:
                    con = max(max(congestion_list + [0]),max([max(E + [0]) for E in B_Edge_congestion]),max(ext_con+[0]))
                    rel_con = con
            ## Obtaining LOC Tree
            LOC = False
            while LOC == False: #Optimización del árbol extendido
                change = False
                switch_count = -1
                for i in range(len(New_BS)):
                    if len(Join_Edges[New_BS[i]])>1: #switch point
                        switch_count = switch_count + 1
                    
                        switched_T = []
                        switched_ext_con = []
                        switched_congestions = []
                        switched_edges = []
                        for j in range(len(Join_Edges[New_BS[i]])): #Para el primer i y j=0 es la primera BFS extensión del árbol (New_T)
                            aux = []
                            switched_ext_con.append(ext_con)
                            switched_edges.append(B_Edges.copy())
                            switched_edges[j][i] = Join_Edges[New_BS[i]][j]
                            switched_T.append(DT + switched_edges[j])
                            switched_congestions.append(B_Edge_congestion.copy()) #No parece necesario copiar todas las congestiones tantas veces
                            switched_DT, VL,VN = spanning_trees.labeling_tree(switched_T[j]) #Esto no es óptimo, este árbol es casi igual al anterior, se debería poder aprovechar
                            if len(L) != 0:
                                for h in range(len(L)):
                                    if L[h][0] == New_BS[i] or L[h][1] == New_BS[i]: #únicas aristas exteriores donde la congestion puede cambiar
                                        v = L[h][0]
                                        w = L[h][1]
                                        tree_path = path_tree_vertices(switched_DT, VL, VN, v, w)
                                        c = sum([edge[2] for edge in tree_path]) + L[h][2]
                                        switched_ext_con[j][h] = c
                            for l in range(len(Join_Edges[New_BS[i]])):
                                if l!= j:
                                    v = Join_Edges[New_BS[i]][l][0]
                                    w = Join_Edges[New_BS[i]][l][1]
                                    tree_path = path_tree_vertices(switched_DT, VL, VN, v, w)
                                    c = sum([edge[2] for edge in tree_path]) + Join_Edges[New_BS[i]][l][2]
                                    aux.append(c)
                            switched_congestions[j][switch_count] = aux # Se cambia la lista de congestiones asociadas a la arista correspondiente (contada por switch_count)
                            if p < np.inf:
                                if len(switched_congestions[j]) == 0:
                                    switch_con = sum([c**p for c in (congestion_list+[0])]) + sum([c**p for c in (switched_ext_con[j]+[0])])
                                    if switch_con != 0:
                                        switch_rel_con = switch_con/(len(congestion_list) + len(ext_con))
                                else:
                                    switch_con = sum([c**p for c in (congestion_list+[0])]) + sum([c**p for c in (switched_ext_con[j]+[0])]) + sum([sum([c**p for c in (E+[0])]) for E in switched_congestions[j]])
                                    switch_rel_con = switch_con/(len(congestion_list) + sum([len(E) for E in B_Edge_congestion])+ len(ext_con))
                            else:
                                switch_con = max(max(congestion_list+[0]),max([max(E+[0]) for E in switched_congestions[j]]), max(switched_ext_con[j]+[0]))
                                switch_rel_con = switch_con
                            if switch_rel_con < rel_con:
                                change = True
                                new_DT = switched_T[j]
                                B_Edge_congestion = switched_congestions[j]
                                B_Edges = switched_edges[j]
                                ext_con = switched_ext_con[j]
                                rel_con = switch_rel_con
                                con = switch_con
                                break
                        if change == True:
                            break
                                
                if change == False:
                    if len(B_Edge_congestion) != 0:
                        for E in B_Edge_congestion:
                            congestion_list= congestion_list + E
                    congestion_list = congestion_list + ext_con
                    DT = new_DT
                    LOC = True  
            if con > last_congestion:
                Check_Out == True
        if Check_Out == False:
            LOC_congestion_lists.append([congestion_list, DT])
            last_congestion = con
    final_cons = []
    for Y in LOC_congestion_lists:
        if p < np.inf:
            final_cons.append(sum([c**p for c in Y[0]]))
        else:
            final_cons.append(max(Y[0]))
    F_C = min(final_cons)
    if p < np.inf:
        FC_C = F_C**(1/p)
    else:
        FC_C = F_C
    index = final_cons.index(F_C)
    print(FC_C, min([max(Y[0]) for Y in LOC_congestion_lists]), min([sum(Y[0]) for Y in LOC_congestion_lists]))
    return FC_C, LOC_congestion_lists[index][0], LOC_congestion_lists[index][1]


def ROC(G,vertex_list, p): #G = planar graph, vertex_list = vertices in plane coordinates (edges assumed to be segments) 
    print('processing graph')
    DG, PL = test_graphs.dual_graph(G,vertex_list) #DG = dual graph, PL = cell list      
    ##
    DG_dict = spanning_trees.graph_as_dictionary(DG)
    ##
    print('end processing')
    last_congestion = np.inf
    last_tree = []
    last_congestion_list = []
    for k in range(len(PL)):
        print('Polygon basis: ', k)
        DT = []
        # ###
        DT_degrees = dict()
        for i in range(len(PL)):
            DT_degrees[i] = [] #for each vertex, we count the weights of incoming edges
        ###
        S = [k]
        CS = [j for j in range(0,k)]+[j for j in range(k+1,len(PL))] 
        congestion_list = []
        L = dict()
        Con_Check = dict()
        for w in DG_dict[k]:
            if not(w[0] in list(L.keys())):
                L[w[0]] = []
            Con_Check[w[0]] = []
            if w[2]==1:
                L[w[0]].append([k,w[0],w[1]])
            else:
                L[w[0]].append([w[0],k,w[1]])
            Con_Check[w[0]] = False #False significa que alguna congestión relativa en la lista L[w[0]] no ha sido calculada
        Check = True
        c_p = 0
        while (len(DT)+1 < len(PL)) and Check:
            WL = []
            CL = []
            Chosen = []
            CKeys = list(L.keys())
            New_Congestions = dict()
            for v in CKeys:
                New_Congestions[v]=[]
            ###### Iniciamos la congestion relativa del subárbol DT
            if p < np.inf and len(congestion_list)!=0:
                rel_con_aux = sum([c**p for c in congestion_list])
            elif p == np.inf and len(congestion_list)!=0:
                rel_con_aux = max(congestion_list)
            elif len(congestion_list)==0:
                rel_con_aux = 0
            #######   
            for vertex in CKeys:
                if len(L[vertex])==1:
                    WL.append(L[vertex][0][2])
                    CL.append(np.inf)
                    Chosen.append(L[vertex][0])
                elif len(L[vertex])>=2: #Computing Relative Congestions
                    WL.append(np.inf)
                    rel_con = rel_con_aux
                    ###############################################
                    Rel_Con = [[] for edge in L[vertex]]
                    rel_con_list = []
                    
                    for l in range(len(L[vertex])):
                        if Con_Check[vertex] == False:
                            Con_Check[vertex] = []
                            new_DT, VL,VN = spanning_trees.labeling_tree(DT+[L[vertex][l]])
                            for j in range(len(L[vertex])):  #Aquí también se puede mejorar la eficiencia de cálculo
                                if j != l: #e != chosen_e:
                                    v = L[vertex][j][0] #e[0]
                                    w = L[vertex][j][1] #e[1]
                                    tree_path = path_tree_vertices(new_DT, VL, VN, v, w)
                                    #print('treepath: ', tree_path)
                                    c = sum([edge[2] for edge in tree_path]) + L[vertex][j][2] #e[2]
                                    Con_Check[vertex].append(c)
                                    Rel_Con[l].append(c)
                                    if p < np.inf:
                                        rel_con = rel_con + c**p
                                    else:
                                        rel_con = max(rel_con, c)
                        else:
                            for c in Con_Check[vertex]:
                                Rel_Con[l].append(c)
                                if p < np.inf:
                                    rel_con = rel_con + c**p
                                else:
                                    rel_con = max(rel_con, c)
                            
                        if p < np.inf:
                            rel_con_list.append(rel_con/(len(congestion_list)+len(Rel_Con[l])))
                        else:
                            rel_con_list.append(rel_con)
                    # print(rel_con_list)
                    ###############################################
                    chosen_index = 0
                    for j in range(1,len(L[vertex])):
                        
                        if rel_con_list[j] < rel_con_list[chosen_index]: #e[2]< chosen_e[2]:
                            
                            chosen_index = j
                        elif rel_con_list[j] == rel_con_list[chosen_index]: #e[2] == chosen_e[2]:        
                            if L[vertex][chosen_index][0]== vertex:
                                vce = L[vertex][chosen_index][1]
                            else:
                                vce = L[vertex][chosen_index][0]                      
                            if L[vertex][j][0]== vertex:
                                ve = L[vertex][j][1]
                            else:
                                ve = L[vertex][j][0]
                            wde = sum(DT_degrees[ve]) #spanning_trees.weighted_degree(DT, ve)
                            wdce = sum(DT_degrees[vce]) #spanning_trees.weighted_degree(DT, vce)
                            
                            if len(DT)==0:
                                de = 1
                                dce = 1
                            else:
                                de = len(DT_degrees[ve]) #spanning_trees.vertex_degree(DT, ve)
                                dce = len(DT_degrees[vce]) #spanning_trees.vertex_degree(DT, vce)
                                
                            if (wde/de) < (wdce/dce):  #Mejorar este control
                                
                                chosen_index = j
                            elif (int(wde) == int(wdce))and (de>dce): #el anterior control ya decide en este caso (revisar)
                                
                                chosen_index = j 
                    
                    New_Congestions[vertex]=Rel_Con[chosen_index]
                    CL.append(rel_con_list[chosen_index])
                    Chosen.append(L[vertex][chosen_index]) 
                else:
                    WL.append(np.inf)
                    CL.append(np.inf)
                    Chosen.append([0,np.inf,1])
            if min(CL) < np.inf:
                chosen_c = CL[0]
                chosen_f = Chosen[0]
                chosen_i = 0
                for i in range(1,len(CL)):
                    if CL[i]< chosen_c or chosen_f[1] == np.inf:
                        chosen_c = CL[i]
                        chosen_f = Chosen[i]
                        chosen_i = i
                    elif CL[i] == chosen_c and Chosen[i][1]!= np.inf:
                        if Chosen[i][0] == CKeys[i]:
                            vce = Chosen[i][1]
                        else:
                            vce = Chosen[i][0]
                        if chosen_f[0] in CS:
                            vfe = chosen_f[1]
                        else:
                            vfe = chosen_f[0]
                        wdfe = sum(DT_degrees[vfe]) #spanning_trees.weighted_degree(DT, vfe)
                        wdce = sum(DT_degrees[vce]) #spanning_trees.weighted_degree(DT, vce)
                        
                        if len(DT)==0:
                            dfe = 1
                            dce = 1
                        else:
                            dfe = len(DT_degrees[vfe]) #spanning_trees.vertex_degree(DT, vfe)
                            dce = len(DT_degrees[vce]) #spanning_trees.vertex_degree(DT, vce)
                            
                        if (wdce/dce) < (wdfe/dfe):
                            chosen_c = CL[i]
                            chosen_f = Chosen[i]
                            chosen_i = i
                        elif (int(wdfe/dfe) == int(wdce/dce))and (dce>dfe):
                            chosen_c = CL[i]
                            chosen_f = Chosen[i]
                            chosen_i = i
            else:
                chosen_c = WL[0]
                chosen_f = Chosen[0]
                chosen_i = 0
                for i in range(1,len(WL)):
                    if WL[i]< chosen_c or chosen_f[1]== np.inf:
                        chosen_c = WL[i]
                        chosen_f = Chosen[i]
                        chosen_i = i
                    elif WL[i] == chosen_c and Chosen[i][1]!= np.inf:
                        if Chosen[i][0] == CKeys[i]:
                            vce = Chosen[i][1]
                        else:
                            vce = Chosen[i][0]
                        if chosen_f[0]in CKeys:
                            vfe = chosen_f[1]
                        else:
                            vfe = chosen_f[0]
                        wdfe = sum(DT_degrees[vfe]) #spanning_trees.weighted_degree(DT, vfe)
                        wdce = sum(DT_degrees[vce]) #spanning_trees.weighted_degree(DT, vce)
                        
                        if len(DT)==0:
                            dfe = 1
                            dce = 1
                        else:
                            dfe = len(DT_degrees[vfe]) #spanning_trees.vertex_degree(DT, vfe)
                            dce = len(DT_degrees[vce]) #spanning_trees.vertex_degree(DT, vce)
                            
                        if (wdce/dce) < (wdfe/dfe):
                            chosen_c = WL[i]
                            chosen_f = Chosen[i]
                            chosen_i = i
                        elif (int(wdfe/dfe) == int(wdce/dce))and (dce>dfe):
                            chosen_c = WL[i]
                            chosen_f = Chosen[i]
                            chosen_i = i
            DT = DT + [chosen_f]
            c_list = New_Congestions[CKeys[chosen_i]]
            congestion_list = congestion_list + c_list
            DT_degrees[chosen_f[0]].append(chosen_f[2])
            DT_degrees[chosen_f[1]].append(chosen_f[2])
            Q = CKeys[chosen_i]
            S.append(Q)
            CS.remove(Q)  #ojo, .remove() puede dar problemas
            #### Updating Congestion #####
            if p<np.inf:
                c_p = c_p + sum([c**p for c in c_list])
            else:
                c_p = max([c_p]+ c_list)
            
            if c_p >= last_congestion:
                Check = False
                
            else:    
                ############ ACTUALIZACION DE L y Con_Check ####################
                L.pop(CKeys[chosen_i]) #ojo, .pop() puede dar problemas
                Con_Check.pop(CKeys[chosen_i])
                for w in DG_dict[Q]:
                    if not(w[0] in S):
                        Con_Check[w[0]] = False #No óptimo si w[0] in CKeys
                        if not(w[0] in CKeys):
                            L[w[0]] = []
                            CKeys.append(w[0])
                        if w[2]==1:
                            L[w[0]].append([Q,w[0],w[1]])
                        else:
                            L[w[0]].append([w[0],Q,w[1]])
                ####################################################
        if (c_p < last_congestion)and(Check == True): 
            last_congestion = c_p
            last_tree = DT
            last_congestion_list = congestion_list
    print(max(last_congestion_list), sum(last_congestion_list), last_congestion)
    print(last_congestion_list)
    return last_congestion, last_congestion_list, last_tree
