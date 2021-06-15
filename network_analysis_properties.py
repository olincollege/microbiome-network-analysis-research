# Network Analysis Properties

import networkx as nx
from networkx.algorithms import community # This part of networkx, for community detection, needs to be imported separately.
from operator import itemgetter

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

def get_network_from_adjm(df_adjm, name=""):
    """
    df_adjm: OTU x OTU adjacency matrix
    """
    G = nx.from_numpy_array(df_adjm.to_numpy())
    G.name = name

    # label nodes with OTUs
    otu_dict = {}
    for node in G.nodes:
        otu_dict[node] = df_adjm.index[node]
    G = nx.relabel_nodes(G, otu_dict)

    # add degree attribute
    degree_dict = dict(G.degree(G.nodes()))
    nx.set_node_attributes(G, degree_dict, 'degree')


    # eigenvalue centrality
    if (name is not "3M"): # power iteration fails to converge within 100 iterations
        eigenvector_dict = nx.eigenvector_centrality(G)
        nx.set_node_attributes(G, eigenvector_dict, 'eigenvector')

    # betweenness centrality
    betweenness_dict = nx.betweenness_centrality(G)
    nx.set_node_attributes(G, betweenness_dict, 'betweenness')

    return G

def get_network_attribute_summary(G):
    """
    G: networkx graph
    """
    # general info
    print(nx.info(G))

    # density
    density = nx.density(G)
    print("Network density:", density)

    # transitivity
    triadic_closure = nx.transitivity(G)
    print("Triadic closure:", triadic_closure)

    # diameter
    # If your Graph has more than one component, this will return False:
    is_connected = nx.is_connected(G)
    print("Network is connected:", is_connected)

    # Next, use nx.connected_components to get the list of components,
    # then use the max() command to find the largest one:
    components = nx.connected_components(G)
    largest_component = max(components, key=len)

    # Create a "subgraph" of just the largest component
    # Then calculate the diameter of the subgraph, just like you did with density.
    subgraph = G.subgraph(largest_component)
    diameter = nx.diameter(subgraph)
    print("Network diameter of largest component:", diameter)


def display_sorted_property(G, attribute="degree"):
    """
    G: networkx graph
    attribute: node attribute
    """
    attr_dict = nx.get_node_attributes(G, attribute)
    sorted_attr = sorted(attr_dict.items(), key=itemgetter(1))
    sorted_attr_OTU = [d[0] for d in sorted_attr] 
    sorted_attr_value = [d[1] for d in sorted_attr]

    f, ax = plt.subplots(figsize=(6, 15))
    ax.barh(sorted_attr_OTU, width=sorted_attr_value)
    ax.set_title(f"{G.name} {attribute}")
    ax.set_xlabel(f"{attribute}")


def visualize_network(G, layout="default"):
    plt.figure()
    if layout == "circular":
        nx.draw_circular(G)
    else:
        nx.draw(G)
    plt.title(f"{G.name} {layout} layout")






