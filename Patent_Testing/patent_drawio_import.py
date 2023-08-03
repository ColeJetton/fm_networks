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


nodes = {}

edges = []

#go line by line
for line in data.splitlines():
   
    #this checks if there's a vertex and adding it to a node
    #only purpose of the node is to ensure that relabeling can occur
    #since the key cell id (in the edge list) is the main thing
    #creates a dictionary of node mapping
    # NOTE: future version should include repeat node testing but
    # idk if that's necessary
    if ('vertex' in line) and ('parent="1"' in line):
        id_index = line.index('id=')
        val_index = line.index('value=')
        style_index = line.index('style=')
        
        nodes.update({line[id_index+4:val_index-2] : line[val_index+7:style_index-2]})

    #essentially 
      
    if ('edge' in line) and not ('edgeLabel' in line) and ('target' in line) and ('source' in line):
        source_index = line.index('source=') 
        initial_index = line[source_index:].index('"')
        
        end_source_index = source_index + 8 + line[source_index+initial_index+1:].index('"')
        
        targ_index = line.index('target=')
        initial_index = line[targ_index:].index('"')
        end_targ_index = targ_index + 8 + line[targ_index+initial_index+1:].index('"')

        #when you come back to this
        edges = edges + [(line[source_index+8:end_source_index], line[targ_index+8:end_targ_index])]

# %%

import networkx as nx

G = nx.Graph()

# I should be able to import it all at once but for some reason it's not working
# so this is the next best thing!
for i in range(0,len(edges)):
    G.add_edge(*edges[i])


#change labeling scheme to prevent repeat name issues
i = 1 
for key in nodes.keys():
    try:
        loc_and = nodes[key].index('&')
        nodes[key] = '(' +  str(i) + ') ' + nodes[key][:loc_and]
    except:
       nodes[key] = '(' +  str(i) + ') ' + nodes[key]
    i+=1
    

G = nx.relabel_nodes(G, nodes)

nx.draw_networkx(G,with_labels=True)


print('EV_C', nx.eigenvector_centrality(G))#,
      #'D_C', nx.degree_centrality(G))    