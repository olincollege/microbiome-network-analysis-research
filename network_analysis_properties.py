# Network Analysis Properties

import networkx as nx
from networkx.algorithms import community # This part of networkx, for community detection, needs to be imported separately.
from operator import itemgetter

import matplotlib.pyplot as plt
from matplotlib import cm

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

    # degree
    degree_dict = dict(G.degree(G.nodes()))
    nx.set_node_attributes(G, degree_dict, 'degree')


    # eigenvalue centrality
    if (name != "3M"): # power iteration fails to converge within 100 iterations
        eigenvector_dict = nx.eigenvector_centrality(G)
        nx.set_node_attributes(G, eigenvector_dict, 'eigenvector_centrality')

    # betweenness centrality
    betweenness_dict = nx.betweenness_centrality(G)
    nx.set_node_attributes(G, betweenness_dict, 'betweenness')

    # calculate modules
    communities = community.greedy_modularity_communities(G)
    M = len(communities)
    modularity_dict = {}
    for i,c in enumerate(communities):
        for name in c:
            modularity_dict[name] = i
    nx.set_node_attributes(G, modularity_dict, 'modularity')

    # within-module degree
    within_module_degree_dict = {}
    for i in G.nodes():
        neighbors = G.neighbors(i)
        module_members = [n for n in G.nodes() if G.nodes[n]['modularity'] == G.nodes[i]['modularity']]
        intersection = [v for v in neighbors if v in module_members]
        within_module_degree_dict[i] = len(intersection)        
    nx.set_node_attributes(G, within_module_degree_dict, 'within_module_degree')

    # within-module connectivity (zi)
    within_module_connectivity_dict = {}
    for i in G.nodes():
        b =  G.nodes[i]['modularity']
        k_ib = G.nodes[i]['within_module_degree']
        bs = [G.nodes[n]['within_module_degree'] for n in G.nodes if G.nodes[n]['modularity'] == b]
        k_bavg = np.mean(bs)
        k_bstd = np.std(bs)
        if (k_bstd == 0):
            within_module_connectivity_dict[i] = 0
        else:
            within_module_connectivity_dict[i] = (k_ib - k_bavg)/(k_bstd)
    nx.set_node_attributes(G, within_module_connectivity_dict, 'within_module_connectivity')

    # among-module connectivity (Pi)
    among_module_connectivity_dict = {}
    for i in G.nodes():
        neighbors = G.neighbors(i)
        k_i = G.nodes[i]['degree']
        sum_k = 0
        for c in range(M):
            module_members = [n for n in G.nodes() if G.nodes[n]['modularity'] == c]
            intersection = [v for v in neighbors if v in module_members]
            k_ic = len(intersection)
            sum_k += (k_ic / k_i) ** 2
        among_module_connectivity_dict[i] = 1 - sum_k
    nx.set_node_attributes(G, among_module_connectivity_dict, 'among_module_connectivity')

    return G


def get_node_roles(G, zt=2.5, Pt=0.62, verbose=False):
    """
    G: networkx graph
    zt: threshold for within-module connectivity
    Pt: threshold for among-module connectivity
    """

    peripheral_nodes = []
    connector_nodes = []
    module_hub_nodes = []
    network_hub_nodes = []

    for i in G.nodes():
        zi = G.nodes[i]['within_module_connectivity']
        Pi = G.nodes[i]['among_module_connectivity']        
        if zi <= zt and Pi <= Pt:
            peripheral_nodes.append((i, zi, Pi))
        elif zi <= zt and Pi > Pt:
            connector_nodes.append((i, zi, Pi))
        elif zi > zt and Pi <= Pt:
            module_hub_nodes.append((i, zi, Pi))
        else:
            network_hub_nodes.append((i, zi, Pi))

    if (verbose):
        print("---")
        print("# peripheral nodes:", len(peripheral_nodes))
        print("# connector nodes:", len(connector_nodes))
        print("# module hub nodes:", len(module_hub_nodes))
        print("# network hub nodes:", len(network_hub_nodes))
        print("total number of nodes:", len(G.nodes))

    node_roles = {
        "peripheral_nodes": peripheral_nodes,
        "connector_nodes": connector_nodes,
        "module_hub_nodes": module_hub_nodes,
        "network_hub_nodes": network_hub_nodes
    }

    return node_roles


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


color_list = [
    "#e377c2",    # pink
    "#f7b6d2",    # light pink
    "#d62728",    # red
    "#ff9896",    # light red
    "#ff7f0e",    # orange
    "#ffbb78",    # light orange
    "#bcbd22",    # mossy green 
    "#dbdb8d",    # light mossy green
    "#8c564b",    # brown
    "#c49c94",    # light brown
    "#2ca02c",    # green
    "#98df8a",    # light green
    "#1f77b4",    # dark blue
    "#aec7e8",    # light blue
    "#9467bd",    # purple
    "#c5b0d5",    # light purple
    "#17becf",    # sky blue
    "#9edae5"     # light sky blue
    "#7f7f7f",    # gray
    "#c7c7c7",    # light gray
]

def visualize_network(G, layout="default"):
    # list of 20 colors
    colors = []
    for i in G.nodes():
        colors.append(color_list[G.nodes[i]["modularity"]])

    plt.figure()
    if layout == "circular":
        nx.draw_circular(G, node_color=colors)
    else:
        nx.draw(G, node_color=colors)
    plt.title(f"{G.name} {layout} layout")


def network_to_csv(G, filename):
    """
    G: networkx graph
    filename: csv to save to
    """
    attributes = set([k for n in G.nodes for k in G.nodes[n].keys()])
    df = pd.DataFrame(index=G.nodes, columns=attributes)
    for i in G.nodes:
        for a in attributes:
            df.loc[i, a] = G.nodes[i][a]
    
    print(f"saving to {filename}")
    df.to_csv(filename)




