# -*- coding: utf-8 -*-
"""
Support functions for functional model loading

Purpose of this file is to include all functions necessary for loading in
previous models and turning them into networks/DSMs that can be analyzed.
@author: colej
"""


def draw_xml_to_usable(xml_file_name):

    # read in 
    with open(xml_file_name, 'r') as f:
        data = f.read()
    
    
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
            print(initial_index)
            print(line[source_index+initial_index+1:].index('"'))
            end_source_index = source_index + 8 + line[source_index+initial_index+1:].index('"')
            print(line)
            
            targ_index = line.index('target=')
            initial_index = line[targ_index:].index('"')
            end_targ_index = targ_index + 8 + line[targ_index+initial_index+1:].index('"')
    
            #when you come back to this
            edges = edges + [(line[source_index+8:end_source_index], line[targ_index+8:end_targ_index])]
    
    
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
        
    
    return G