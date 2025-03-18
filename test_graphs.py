import random as rd
import numpy as np
import spanning_trees
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import collections  as mc

def complete_graph(n):
    graph = []
    for i in range(n):
        if i < (n-1):
            for j in range(i+1,n):
                graph.append([i,j,1])
    return graph

def complete_graph_sum_weighted(n):
    graph = []
    for i in range(n):
        if i < (n-1):
            for j in range(i+1,n):
                graph.append([i,j,i+j])
    return graph

def complete_graph_dif_weighted(n):
    graph = []
    for i in range(n):
        if i < (n-1):
            for j in range(i+1,n):
                graph.append([i,j,abs(j-i)])
    return graph

def complete_bipartite_graph(n,m):
    graph = []
    for i in range(n):
        for j in range(n,m+n):
            graph.append([i,j,1])
    return graph

def complete_bipartite_graph_sum_weighted(n,m):
    graph = []
    for i in range(n):
        for j in range(n,m+n):
            graph.append([i,j,i+j])
    return graph

def complete_bipartite_graph_dif_weighted(n,m):
    graph = []
    for i in range(n):
        for j in range(n,m+n):
            graph.append([i,j,abs(j-i)])
    return graph


def complete_multipartite(order_list):
    number_sets = len(order_list)
    graph = []
    for i in range(number_sets-1):
        actual_index = sum(order_list[0:i])
        for j in range(i+1,number_sets):
            if j != i:
                moving_index = sum(order_list[0:j])
                for k in range(order_list[i]):
                    for l in range(order_list[j]):
                        graph.append([actual_index + k , moving_index + l,1])
    return graph

def complete_multipartite_sum_weighted(order_list):
    number_sets = len(order_list)
    graph = []
    for i in range(number_sets-1):
        actual_index = sum(order_list[0:i])
        for j in range(i+1,number_sets):
            if j != i:
                moving_index = sum(order_list[0:j])
                for k in range(order_list[i]):
                    for l in range(order_list[j]):
                        graph.append([actual_index + k , moving_index + l, actual_index + k + moving_index + l])
    return graph

def complete_multipartite_dif_weighted(order_list):
    number_sets = len(order_list)
    graph = []
    for i in range(number_sets-1):
        actual_index = sum(order_list[0:i])
        for j in range(i+1,number_sets):
            if j != i:
                moving_index = sum(order_list[0:j])
                for k in range(order_list[i]):
                    for l in range(order_list[j]):
                        graph.append([actual_index + k , moving_index + l, abs(actual_index + k - (moving_index + l))])
    return graph

def random_graph(number_vertices, edge_umbral):
    RG = []
    for i in range(number_vertices):
        for j in range(i):
            p = rd.random()
            if p < edge_umbral:
                #random_weight = rd.choice(range(1,10)) #weights between 1 and 9
                RG.append([i,j,1]) #random_weight])
    return RG

def random_graph_sum_weighted(number_vertices, edge_umbral):
    RG = []
    for i in range(number_vertices):
        for j in range(i):
            p = rd.random()
            if p < edge_umbral:
                #random_weight = rd.choice(range(1,10)) #weights between 1 and 9
                RG.append([i,j,i+j]) #random_weight])
    return RG

def grid_graph(n,m):
    GG = []
    vertex_list = []
    for j in range(m):
        for i in range(n):
            vertex_list.append(np.array([i,j]))
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,1])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),1])
    return GG, vertex_list

def grid_graph_sum_weighted(n,m):
    GG = []
    vertex_list = []
    for j in range(m):
        for i in range(n):
            vertex_list.append(np.array([i,j]))
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,n*j + i + n*(j+1)+i])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),n*j+i + (i+1+n*j)])
    return GG, vertex_list

def grid_graph_dif_weighted(n,m):
    GG = []
    vertex_list = []
    for j in range(m):
        for i in range(n):
            vertex_list.append(np.array([i,j]))
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,abs(n*j + i - (n*(j+1)+i))])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),abs(n*j+i - (i+1+n*j))])
    return GG, vertex_list

def cyclic_grid(n,m):
    GG = []
    for i in range(n):
        for j in range(m):
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,1])
            elif j == m-1:
                GG.append([n*j + i,i,1])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),1])
            elif i == n-1:
                GG.append([n*j+i, n*j,1])
    return GG

def grid3D_graph(n,m,k):
    GG = []
    for i in range(n):
        for j in range(m):
            for l in range(k):
                if l != k-1:
                    GG.append([l*(n*m) + n*j + i,(l+1)*(n*m) + n*j + i,1])
                if j != m-1:
                    GG.append([l*(n*m)+n*j + i,l*(n*m)+n*(j+1)+i,1])
                if i != n-1:
                    GG.append([l*(n*m)+n*j+i,l*(n*m)+(i+1+n*j),1])
    return GG

def triangular_grid_graph(k):
    GG = []
    vl =[]
    row_index = 0
    for i in range(k-1):
        for j in range(k-i-1):
           
            if not([j,i]in vl):
                vl.append([j,i])
            GG.append([j + row_index ,j+ row_index+1,1])
            GG.append([j+row_index, j+row_index + (k-i),1])
            GG.append([j+1 + row_index, j + row_index + (k-i),1])
            if not([j+1,i]in vl):
                vl.append([j+1,i])
        row_index = row_index + (k-i)
    vl.append([0,k-1])
    vl = [np.array(v) for v in vl]
    return GG,  vl

def triangular_grid_graph_sum_weighted(k):
    GG = []
    vl =[]
    row_index = 0
    for i in range(k-1):
        for j in range(k-i-1):
           
            if not([j,i]in vl):
                vl.append([j,i])
            GG.append([j + row_index ,j+ row_index+1,j + row_index +j+ row_index+1])
            GG.append([j+row_index, j+row_index + (k-i),j+row_index+ j+row_index + (k-i)])
            GG.append([j+1 + row_index, j + row_index + (k-i),j+1 + row_index + j + row_index + (k-i)])
            if not([j+1,i]in vl):
                vl.append([j+1,i])
        row_index = row_index + (k-i)
    vl.append([0,k-1])
    vl = [np.array(v) for v in vl]
    return GG,  vl

def triangular_grid_graph_dif_weighted(k):
    GG = []
    vl =[]
    row_index = 0
    for i in range(k-1):
        for j in range(k-i-1):
           
            if not([j,i]in vl):
                vl.append([j,i])
            GG.append([j + row_index ,j+ row_index+1,abs(j + row_index -(j+ row_index+1))])
            GG.append([j+row_index, j+row_index + (k-i),abs(j+row_index - (j+row_index + (k-i)))])
            GG.append([j+1 + row_index, j + row_index + (k-i),abs(j+1 + row_index - (j + row_index + (k-i)))])
            if not([j+1,i]in vl):
                vl.append([j+1,i])
        row_index = row_index + (k-i)
    vl.append([0,k-1])
    vl = [np.array(v) for v in vl]
    return GG,  vl

### Planar Graphs ###
            
def random_planar_graph(n): #generates a planar geodesic graph of n random vertices in [0,1]x[0,1]
    GG=[]
    GC = []
    vertex_list = []
    for i in range(n):
        vertex = np.array([rd.random(),rd.random()])
        #print(vertex)
        m = len(GG)
        k = len(vertex_list)
        vertex_list.append(vertex)
        if k != 0:
            if m != 0:
                for l in range(k):
                    vertex_c = vertex_list[l]
                    token = True
                    count = 0
                    while (token == True)&(count < m):
                        a = vertex_list[GG[count][0]]
                        b = vertex_list[GG[count][1]]
                        if (GG[count][0] != l) and (GG[count][1] != l):
                            u = vertex-a
                            v = vertex_c-a
                            w = b-a
                            det_uw = u[0]*w[1] - u[1]*w[0]
                            det_vw = v[0]*w[1] - v[1]*w[0]
                            # dot_uw = np.dot(u,w)
                            # dot_vw = np.dot(v,w)
                            #print(det_uw,det_vw)
                            x = a-vertex 
                            y = b-vertex
                            z = vertex_c - vertex
                            det_xz = x[0]*z[1] - x[1]*z[0]
                            det_yz = y[0]*z[1] - y[1]*z[0]
                            # dot_xz = np.dot(x,z)
                            # dot_yz = np.dot(y,z)
                            #print(det_xz,det_yz)
                            if (det_uw*det_vw < 0) and (det_xz*det_yz < 0): #and (dot_uw < 0) and (dot_vw < 0)  and (dot_xz < 0) and (dot_yz < 0):
                                token = False                                
                        count = count + 1
                    if token == True:
                        GG.append([l,k,1])# np.linalg.norm(vertex-vertex_c)])
                        GC.append([l,k])
            else:
                GC.append([0,1])
                GG.append([0,1,1]) #np.linalg.norm(vertex-vertex_list[0])])
    
    lines = [[vertex_list[e[0]],vertex_list[e[1]]] for e in GG]
    lc = mc.LineCollection(lines, linewidths=2)
    fig, ax = pl.subplots()
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
    plt.show()
    return GG, vertex_list   

def random_planar_graph_metric_weighted(n): #generates a planar geodesic graph of n random vertices in [0,1]x[0,1]
    GG=[]
    GC = []
    vertex_list = []
    for i in range(n):
        vertex = np.array([rd.random(),rd.random()])
        #print(vertex)
        m = len(GG)
        k = len(vertex_list)
        vertex_list.append(vertex)
        if k != 0:
            if m != 0:
                for l in range(k):
                    vertex_c = vertex_list[l]
                    token = True
                    count = 0
                    while (token == True)&(count < m):
                        a = vertex_list[GG[count][0]]
                        b = vertex_list[GG[count][1]]
                        if (GG[count][0] != l) and (GG[count][1] != l):
                            u = vertex-a
                            v = vertex_c-a
                            w = b-a
                            det_uw = u[0]*w[1] - u[1]*w[0]
                            det_vw = v[0]*w[1] - v[1]*w[0]
                            # dot_uw = np.dot(u,w)
                            # dot_vw = np.dot(v,w)
                            #print(det_uw,det_vw)
                            x = a-vertex 
                            y = b-vertex
                            z = vertex_c - vertex
                            det_xz = x[0]*z[1] - x[1]*z[0]
                            det_yz = y[0]*z[1] - y[1]*z[0]
                            # dot_xz = np.dot(x,z)
                            # dot_yz = np.dot(y,z)
                            #print(det_xz,det_yz)
                            if (det_uw*det_vw < 0) and (det_xz*det_yz < 0): #and (dot_uw < 0) and (dot_vw < 0)  and (dot_xz < 0) and (dot_yz < 0):
                                token = False                                
                        count = count + 1
                    if token == True:
                        GG.append([l,k,np.linalg.norm(vertex-vertex_c)]) # 1])# 
                        GC.append([l,k])
            else:
                GC.append([0,1])
                GG.append([0,1,np.linalg.norm(vertex-vertex_list[0])]) # 1]) #
    
    lines = [[vertex_list[e[0]],vertex_list[e[1]]] for e in GG]
    lc = mc.LineCollection(lines, linewidths=2)
    fig, ax = pl.subplots()
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
    plt.show()
    return GG, vertex_list  

def double_random_planar_graph(n,p): #generates a planar geodesic graph of n random vertices in [0,1]x[0,1], each new edge with probability p of appearance
    GG=[]
    vertex_list = []
    for i in range(n):
        vertex = np.array([rd.random(),rd.random()])
        #print(vertex)
        m = len(GG)
        k = len(vertex_list)
        vertex_list.append(vertex)
        if k != 0:
            if m != 0:
                for l in range(k):
                    vertex_c = vertex_list[l]
                    token = True
                    count = 0
                    while (token == True)&(count < m):
                        a = vertex_list[GG[count][0]]
                        b = vertex_list[GG[count][1]]
                        if (GG[count][0] != l) and (GG[count][1] != l):
                            u = vertex-a
                            v = vertex_c-a
                            w = b-a
                            det_uw = u[0]*w[1] - u[1]*w[0]
                            det_vw = v[0]*w[1] - v[1]*w[0]
                            x = a-vertex 
                            y = b-vertex
                            z = vertex_c - vertex
                            det_xz = x[0]*z[1] - x[1]*z[0]
                            det_yz = y[0]*z[1] - y[1]*z[0]
                            if (det_uw*det_vw < 0) and (det_xz*det_yz < 0): #and (dot_uw < 0) and (dot_vw < 0)  and (dot_xz < 0) and (dot_yz < 0):
                                token = False                                
                        count = count + 1
                    if token == True:
                        q = rd.random()
                        if q < p:
                            GG.append([l,k,np.linalg.norm(vertex-vertex_c)])
            else:
                q = rd.random()
                if q < p:
                    GG.append([0,1,np.linalg.norm(vertex-vertex_list[0])])
    return GG, vertex_list  

def quad_triangle_grid(n,m):
    GG = []
    vertex_list = []
    for j in range(m):
        for i in range(n):
            vertex_list.append(np.array([i,j]))
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,1])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),1])
    for j in range(m-1):
        for i in range(n-1):
            if (j)%2 == 0:
                GG.append([ n*j+i, n*(j+1)+i+1,1])
            else:
                GG.append([ n*(j+1)+i, n*(j)+i+1,1])
    return GG, vertex_list

def quad_triangle_grid_sum_weighted(n,m):
    GG = []
    vertex_list = []
    for j in range(m):
        for i in range(n):
            vertex_list.append(np.array([i,j]))
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,n*j + i + n*(j+1)+i])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),n*j+i + (i+1+n*j)])
    for j in range(m-1):
        for i in range(n-1):
            if (j)%2 == 0:
                GG.append([ n*j+i, n*(j+1)+i+1, n*j+i + n*(j+1)+i+1])
            else:
                GG.append([ n*(j+1)+i, n*(j)+i+1, n*(j+1)+i, n*(j)+i+1])
    return GG, vertex_list

def quad_triangle_grid_diagonals(n,m):
    GG = []
    vertex_list = []
    for j in range(m):
        for i in range(n):
            vertex_list.append(np.array([i,j]))
            if j != m-1:
                GG.append([n*j + i,n*(j+1)+i,1])
            if i != n-1:
                GG.append([n*j+i,(i+1+n*j),1])
    for j in range(m-1):
        for i in range(n-1):
            GG.append([ n*j+i, n*(j+1)+i+1,1])
    return GG, vertex_list

def triangle_hex_grid(n):
    K = 2*n+1
    GG = []
    vertex_list = []
    for j in range(n+1):
        if j==0:
            for i in range(0,K):       
                vertex_list.append(np.array([i ,j]))
                if  j < n:
                    if (K*j + i)%2 == 0: # and j%2 == 0:
                        GG.append([K*j + i,K*(j+1)+i,1])
                    # elif (K*j + i)%2 == 0 and j%2 ==1:
                    #     GG.append([K*j + i-int(j/2),K*(j+1)+i-int((j+1)/2),1])
                    if i != K-1:
                        GG.append([K*j+i,(i+1+K*j),1])
                elif j==n:
                    if i != K-1:
                        GG.append([K*j+i,(i+1+K*j),1])
        else:
            for i in range(j-1,K-j+1):       #corregir índice, sigue un número triangular
                vertex_list.append(np.array([i ,j]))
                if  j < n:
                    if (K*j + i)%2 == 0: # and j%2 == 0:
                        GG.append([K*j + i-2*int(((j-1)*(j-2))/2)-j+1,K*(j+1)+i-2*int(((j)*(j-1))/2)-j,1])
                    # elif (K*j + i)%2 == 0 and j%2 ==1:
                    #     GG.append([K*j + i-int(j/2),K*(j+1)+i-int((j+1)/2),1])
                    if i != K-j:
                        GG.append([K*j+i-2*int(((j-1)*(j-2))/2)-j+1,(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1,1])
                elif j==n:
                    if i != K-j:
                        GG.append([K*j+i-2*int(((j-1)*(j-2))/2)-j+1,(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1,1])
    return GG, vertex_list

def triangle_hex_grid_sum_weighted(n):
    K = 2*n+1
    GG = []
    vertex_list = []
    for j in range(n+1):
        if j==0:
            for i in range(0,K):       
                vertex_list.append(np.array([i ,j]))
                if  j < n:
                    if (K*j + i)%2 == 0: # and j%2 == 0:
                        GG.append([K*j + i,K*(j+1)+i,K*j + i+K*(j+1)+i])
                    # elif (K*j + i)%2 == 0 and j%2 ==1:
                    #     GG.append([K*j + i-int(j/2),K*(j+1)+i-int((j+1)/2),1])
                    if i != K-1:
                        GG.append([K*j+i,(i+1+K*j),K*j+i+(i+1+K*j)])
                elif j==n:
                    if i != K-1:
                        GG.append([K*j+i,(i+1+K*j),K*j+i+(i+1+K*j)])
        else:
            for i in range(j-1,K-j+1):       #corregir índice, sigue un número triangular
                vertex_list.append(np.array([i ,j]))
                if  j < n:
                    if (K*j + i)%2 == 0: # and j%2 == 0:
                        GG.append([K*j + i-2*int(((j-1)*(j-2))/2)-j+1,K*(j+1)+i-2*int(((j)*(j-1))/2)-j,K*j + i-2*int(((j-1)*(j-2))/2)-j+1+K*(j+1)+i-2*int(((j)*(j-1))/2)-j])
                    # elif (K*j + i)%2 == 0 and j%2 ==1:
                    #     GG.append([K*j + i-int(j/2),K*(j+1)+i-int((j+1)/2),1])
                    if i != K-j:
                        GG.append([K*j+i-2*int(((j-1)*(j-2))/2)-j+1,(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1,K*j+i-2*int(((j-1)*(j-2))/2)-j+1+(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1])
                elif j==n:
                    if i != K-j:
                        GG.append([K*j+i-2*int(((j-1)*(j-2))/2)-j+1,(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1,K*j+i-2*int(((j-1)*(j-2))/2)-j+1+(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1])
    return GG, vertex_list

def triangle_hex_grid_dif_weighted(n):
    K = 2*n+1
    GG = []
    vertex_list = []
    for j in range(n+1):
        if j==0:
            for i in range(0,K):       
                vertex_list.append(np.array([i ,j]))
                if  j < n:
                    if (K*j + i)%2 == 0: # and j%2 == 0:
                        GG.append([K*j + i,K*(j+1)+i,abs(K*j + i-(K*(j+1)+i))])
                    # elif (K*j + i)%2 == 0 and j%2 ==1:
                    #     GG.append([K*j + i-int(j/2),K*(j+1)+i-int((j+1)/2),1])
                    if i != K-1:
                        GG.append([K*j+i,(i+1+K*j),abs(K*j+i - (i+1+K*j))])
                elif j==n:
                    if i != K-1:
                        GG.append([K*j+i,(i+1+K*j),abs(K*j+i - (i+1+K*j))])
        else:
            for i in range(j-1,K-j+1):       #corregir índice, sigue un número triangular
                vertex_list.append(np.array([i ,j]))
                if  j < n:
                    if (K*j + i)%2 == 0: # and j%2 == 0:
                        GG.append([K*j + i-2*int(((j-1)*(j-2))/2)-j+1,K*(j+1)+i-2*int(((j)*(j-1))/2)-j,abs(K*j + i-2*int(((j-1)*(j-2))/2)-j+1 - (K*(j+1)+i-2*int(((j)*(j-1))/2)-j))])
                    # elif (K*j + i)%2 == 0 and j%2 ==1:
                    #     GG.append([K*j + i-int(j/2),K*(j+1)+i-int((j+1)/2),1])
                    if i != K-j:
                        GG.append([K*j+i-2*int(((j-1)*(j-2))/2)-j+1,(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1,abs(K*j+i-2*int(((j-1)*(j-2))/2)-j+1 - ((i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1))])
                elif j==n:
                    if i != K-j:
                        GG.append([K*j+i-2*int(((j-1)*(j-2))/2)-j+1,(i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1,abs(K*j+i-2*int(((j-1)*(j-2))/2)-j+1 - ((i+1+K*j)-2*int(((j-1)*(j-2))/2)-j+1))])
    return GG, vertex_list

def quad_hex_grid(n,m):
    k = 2*(n+1)
    l = m + 1
    GG = []
    vertex_list = []
    for j in range(l):
        for i in range(k):
            
            if 0 < j and j < l-1:
                vertex_list.append(np.array([i,j]))
                if j != l-2 and (k*j + i)%2 == 0 and j%2 == 0:
                    GG.append([k*j + i-1,k*(j+1)+i-1,1])
                elif j != l-2 and (k*j + i)%2 == 1 and j%2 ==1:
                    GG.append([k*j + i-1,k*(j+1)+i-1,1])
                if j == l-2 and (k*j + i)%2 == 0 and j%2 == 0:
                    GG.append([k*j + i-1,k*(j+1)+i-1,1])
                elif j == l-2 and (k*j + i)%2 == 1 and j%2 ==1:
                    GG.append([k*j + i-1,k*(j+1)+i-2,1])
                if i != k-1:
                    GG.append([k*j+i-1,(i+k*j),1])
            elif j == 0 and i!= k-1:
                vertex_list.append(np.array([i,j]))
                if j != l-1 and (k*j + i)%2 == 0 and j%2 == 0:
                    GG.append([k*j + i,k*(j+1)+i-1,1])
                elif j != l-1 and (k*j + i)%2 == 1 and j%2 ==1:
                    GG.append([k*j + i,k*(j+1)+i-1,1])
                if i != k-2:
                    GG.append([k*j+i,(i+1+k*j),1])
            elif (l%2 == 1 and j == l-1 and i!= 0):
                vertex_list.append(np.array([i,j]))
                if i != k-1:
                    GG.append([k*j+i-2,(i-1+k*j),1])
            elif (l%2 == 0 and j == l-1 and i!= k-1):
                vertex_list.append(np.array([i,j]))
                if i != k-2:
                    GG.append([k*j+i-1,(i+k*j),1])
            
    return GG, vertex_list

#####################################################################

def view_planar_graph(Planar_Graph, Vertex_List):
    lines = [[Vertex_List[e[0]],Vertex_List[e[1]]] for e in Planar_Graph]
    lc = mc.LineCollection(lines, linewidths=2)
    
    fig, ax = pl.subplots()
    for v in Vertex_List:
        plt.scatter(v[0],v[1], color = 'black', zorder = 2)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
    plt.show()
    return None

def hypercube_graph(n): #n is the dimension of the hypercube
    G = [[0,1,1]]
    for i in range(1,n):
        G_copy = G.copy()
        for e in G_copy:
            G.append([e[0]+2**i,e[1]+2**i,1])
        for j in range(2**i,2**(i+1)):
                G.append([j,j-2**i,1])
    return G

def cycle_ordering(Planar_Graph, Vertex_List, Vertex_Index):
    star = spanning_trees.vertex_star_vertex(Planar_Graph, Vertex_Index)
    arguments = []
    for i in star:
        aux = Vertex_List[i] - Vertex_List[Vertex_Index]
        complex_number = aux[0] + aux[1]*1j
        arguments.append(np.angle(complex_number))
    sorted_indices = np.argsort(arguments)
    ordered_list = []
    for j in range(len(star)):
        ordered_list.append(star[sorted_indices[j]])    
    return ordered_list

def consecutive_clockwise_vertex(Planar_Graph,Vertex_List, i, j): #[i,j] or [j,i] is an edge
    star = cycle_ordering(Planar_Graph, Vertex_List, i)
    index = star.index(j)
    cons_vertex = star[index - 1]
    return cons_vertex

def consecutive_counterclockwise_vertex(Planar_Graph,Vertex_List, i, j): #[i,j] or [j,i] is an edge
    star = cycle_ordering(Planar_Graph, Vertex_List, i)
    index = star.index(j)
    m = len(star)
    cons_vertex = star[(index + 1)%m]
    return cons_vertex

def polygons(Planar_Graph,Vertex_List, edge): #Planar graph must be reduced
    polygon_1 = {edge[0],edge[1]}
    polygon_2 = {edge[0],edge[1]}
    counter_token = True
    u = edge[0]
    w = edge[1]
    while counter_token == True:
        v = consecutive_counterclockwise_vertex(Planar_Graph, Vertex_List, w, u)
        if not(v in polygon_1):
            polygon_1.add(v)
            u = w
            w = v
        else:
            counter_token = False
    clock_token = True
    u = edge[0]
    w = edge[1]
    while clock_token == True:
        v = consecutive_clockwise_vertex(Planar_Graph, Vertex_List, w, u)
        if not(v in polygon_2):
            polygon_2.add(v)
            u = w
            w = v
        else:
            clock_token = False
    return polygon_1, polygon_2

def dual_graph(Planar_Graph, Vertex_List): #bastante ineficiente pero funciona #pueden aparecer aristas dobles si el grafo no está reducido
    GG = []
    Dual_Polygons = []
    count = 0
    for edge in Planar_Graph:
        P1,P2 = polygons(Planar_Graph, Vertex_List, edge)
        if P1 in Dual_Polygons:
            index_P1 = Dual_Polygons.index(P1)
        else:
            Dual_Polygons.append(P1)
            index_P1 = count
            count = count + 1
        if P2 in Dual_Polygons:
            index_P2 = Dual_Polygons.index(P2)
        else:
            Dual_Polygons.append(P2)
            index_P2 = count
            count = count + 1
        GG.append([index_P1,index_P2,edge[2]])
    return GG, Dual_Polygons

def Kruskal_min(GG):
    Kruskal = []
    Components = []
    Vertices = []
    GS = sorted(GG, key=lambda tup: tup[2])
    for e in GS:
        if not(e[0] in Vertices) and not(e[1] in Vertices):
            Vertices.append(e[0])
            Vertices.append(e[1])
            Kruskal.append(e)
            Components.append({e[0],e[1]})
        
        elif e[0] in Vertices and not(e[1] in Vertices):
            Vertices.append(e[1])
            Kruskal.append(e)
            for C in Components:
                if e[0] in C:
                    C.add(e[1])
        elif e[1] in Vertices and not(e[0] in Vertices):
            Vertices.append(e[0])
            Kruskal.append(e)
            for C in Components:
                if e[1] in C:
                    C.add(e[0])
        elif e[0] in Vertices and e[1] in Vertices:
            for C in Components:
                if e[0] in C:
                    index_0 = Components.index(C)
                if e[1] in C:
                    index_1 = Components.index(C)
            if index_0 != index_1:
                Kruskal.append(e)
                Components[index_0] = Components[index_0].union(Components[index_1])
                Components.remove(Components[index_1])
    return Kruskal

def first_search_rooted_tree(GG,v):
    FS_Tree = []
    Vertices = {v}
    Boundary_Vertices = {v}
    N = spanning_trees.cardinal(GG)
    while len(Vertices)<N:
        new_Boundary_Vertices = set()
        for u in Boundary_Vertices:
            S = spanning_trees.vertex_star(GG,u)
            for e in S:
                if not(e[0] in Vertices):
                    Vertices.add(e[0])
                    new_Boundary_Vertices.add(e[0])
                    FS_Tree.append(e)
                elif not(e[1] in Vertices):
                    Vertices.add(e[1])
                    new_Boundary_Vertices.add(e[1])
                    FS_Tree.append(e)
        Boundary_Vertices = new_Boundary_Vertices
    return FS_Tree           

def dual_tree(Planar_Graph, Vertex_List, Dual_Graph, Polygon_List, dual_spanning_tree): #Grafo plano reducido
    forbidden_edges = []
    for edge in dual_spanning_tree:
        edge_intersection = Polygon_List[edge[0]].intersection(Polygon_List[edge[1]])
        if len(edge_intersection) == 2:
            forbidden_edges.append(edge_intersection)
    dual_tree = []
    for edge in Planar_Graph:
        if not({edge[0],edge[1]} in forbidden_edges):
            dual_tree.append(edge)
    dual_tree = Kruskal_min(dual_tree)
    return dual_tree