# -*- coding: utf-8 -*-
"""
@author: colej

requires:
    beautifulsoup4 for text parsing
"""

import pandas as pd


# Reading data from the xml file
patent_name = 'Patent 2982 BTA v1.xml'
#patent_name = 'Patent 7776 BTA v5.xml'
#patent_name = 'Functional Model 7776 v2.xml'
#patent_name = 'Patent 10975 BTA v3.xml'

with open(patent_name, 'r') as f:
    data = f.read()
    

"""

IF it says “vertex”, it is a node.
If it is not, then it is a connection edge (also it says edge). 
They will all include a source and a target. 


text = bs_data.prettify()

"""

# %%
# creating the edgelist!!!!


fxn_nodes = {}
flw_nodes = {}
edges = []

#go line by line
for line in data.splitlines():
   
    #this checks if there's a vertex and adding it to a node
    #only purpose of the node is to ensure that relabeling can occur
    #since the key cell id (in the edge list) is the main thing
    #creates a dictionary of node mapping
    # NOTE: future version should include repeat node testing but
    # idk if that's necessary
    if 'vertex' in line:
        id_index = line.index('id=')
        val_index = line.index('value=')
        style_index = line.index('style=')
        
        fxn_nodes.update({line[id_index+4:val_index-2] : line[val_index+7:style_index-2]})

    #essentially need to record each edge pair then use that with the key later
    #to connect and then relabel them
    
    #however, it looks as if I can create a bipartite graph using the edges 
    #and the names of the edges...
    
    """
    psuedo
    if edge in line
        do exact same thing as before,
        if theres a target
            use indexing to connect to target
        if there's a source
            use indexing to connect to source
    
    
    """
    
    
    # Original (not sure why I included "and not edgelabel" makes no sense)
    # if ('edge' in line) and not ('edgeLabel' in line) and ('target' in line) and ('source' in line):
    
    if 'edge' in line:

        id_index = line.index('id=')
        val_index = line.index('value=')
        style_index = line.index('style=')
        
        # add dictionary for mapping
        flw_nodes.update({line[id_index+4:val_index-2] : line[val_index+7:style_index-2]})
        
        
        
        if 'source' in line:
            source_index = line.index('source=')
            initial_index = line[source_index:].index('"')            
            end_source_index = source_index + 8 + line[source_index+initial_index+1:].index('"')
            #add source to flow edge
            edges = edges + [(line[source_index+8:end_source_index], line[id_index+4:val_index-2] )]
            
        if 'target' in line:
            targ_index = line.index('target=')
            initial_index = line[targ_index:].index('"')
            end_targ_index = targ_index + 8 + line[targ_index+initial_index+1:].index('"')
            edges = edges + [(line[targ_index+8:end_targ_index], line[id_index+4:val_index-2] )]

# %%

import networkx as nx

G = nx.Graph()

# I should be able to import it all at once but for some reason it's not working
# so this is the next best thing!
for i in range(0,len(edges)):
    G.add_edge(*edges[i])


#change labeling scheme to prevent repeat name issues



i = 1 
for key in fxn_nodes.keys():
    try:
        loc_and = fxn_nodes[key].index('&')
        fxn_nodes[key] = '(' +  str(i) + ') ' + fxn_nodes[key][:loc_and]
    except:
       fxn_nodes[key] = '(' +  str(i) + ') ' + fxn_nodes[key]
    i+=1

    
for key in flw_nodes.keys():
    try:
        loc_and = flw_nodes[key].index('&')
        flw_nodes[key] = '(' +  str(i) + ') ' + flw_nodes[key][:loc_and]
    except:
       flw_nodes[key] = '(' +  str(i) + ') ' + flw_nodes[key]
    i+=1


final_dict = {**fxn_nodes, **flw_nodes}
"""
color_map = []
for key in final_dict:
    print(key)
    if key in fxn_nodes:
        print(key)
        color_map.append('blue')
    else:
        color_map.append('green')
"""
G = nx.relabel_nodes(G, final_dict)

color_map = []
for node in G:
    if node in fxn_nodes.values():
        color_map.append('#f7ba63') #light orange for fxns
    else:
        color_map.append('#96b7fa') #light blue for flws
        


import matplotlib.pyplot as plt

#positions = nx.spring_layout(G) #default layout, comment out for others

#positions = nx.planar_layout(G) #non intersecting, looks a bit odd though
positions =nx.kamada_kawai_layout(G) #looks the most like the model! Use this

plt.figure()   
nx.draw_networkx(G , pos = positions, node_color = color_map, with_labels=True, font_size = 8)


#print('EV_C', nx.eigenvector_centrality(G))#,
      #'D_C', nx.degree_centrality(G))    
      
Q = nx.adjacency_matrix(G)
Q = Q.todense()      

plt.figure()
plt.imshow(Q)