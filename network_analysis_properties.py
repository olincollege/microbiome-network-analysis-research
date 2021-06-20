# Network Analysis Properties

import networkx as nx
import networkx.algorithms.community as nx_comm # This part of networkx, for community detection, needs to be imported separately.
from operator import itemgetter

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib import ticker

import pandas as pd
import numpy as np

from collections import Counter

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
    if (name != "3M"):
        eigenvector_dict = nx.eigenvector_centrality(G)
    else:
        # power iteration fails to converge within 100 iterations
        # use alternate implementation
        eigenvector_dict = nx.eigenvector_centrality_numpy(G)
    nx.set_node_attributes(G, eigenvector_dict, 'eigenvector_centrality')

    # betweenness centrality
    betweenness_dict = nx.betweenness_centrality(G)
    nx.set_node_attributes(G, betweenness_dict, 'betweenness')

    # calculate modules
    communities = nx_comm.greedy_modularity_communities(G)
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

    # node roles
    zt=2.5
    Pt=0.62
    node_roles_dict = {}
    for i in G.nodes():
        zi = G.nodes[i]['within_module_connectivity']
        Pi = G.nodes[i]['among_module_connectivity']
        if zi <= zt and Pi <= Pt:
            node_roles_dict[i] = "peripheral_node"
        elif zi <= zt and Pi > Pt:
            node_roles_dict[i] = "connector_node"
        elif zi > zt and Pi <= Pt:
            node_roles_dict[i] = "module_hub_node"
        else:
            node_roles_dict[i] = "network_hub_node"
    nx.set_node_attributes(G, node_roles_dict, 'node_role')

    return G


def networks_attributes_to_csv(Gs, filename):
    """
    Gs: list of Graphs
    """
    names = [G.name for G in Gs]
    df = pd.DataFrame(index=names)

    for G in Gs:
        df.loc[G.name, "number_of_nodes"] = len(G.nodes)

        df.loc[G.name, "number_of_edges"] = len(G.edges)

        df.loc[G.name, "density"] = nx.density(G)

        df.loc[G.name, "transitivity"] = nx.transitivity(G)

        df.loc[G.name, "network_is_connected"] = nx.is_connected(G)

        # network diameter of largest component
        # use nx.connected_components to get the list of components
        components = nx.connected_components(G)
        largest_component = max(components, key=len)
        subgraph = G.subgraph(largest_component)
        # Then calculate the diameter of the subgraph, just like you did with density.
        df.loc[G.name, "network_diameter_of_largest_component"] = nx.diameter(subgraph)
        
        # modularity
        communities = nx_comm.greedy_modularity_communities(G)
        df.loc[G.name, "modularity_score"] = nx_comm.modularity(G, communities)
        df.loc[G.name, "number_of_modules"] = len(communities)

        # number of node roles
        node_role_counts = Counter(dict(G.nodes(data="node_role")).values())
        df.loc[G.name, "number_of_peripheral_node"] = node_role_counts["peripheral_node"]
        df.loc[G.name, "number_of_connector_node"] = node_role_counts["connector_node"]
        df.loc[G.name, "number_of_module_hub_node"] = node_role_counts["module_hub_node"]
        df.loc[G.name, "number_of_network_hub_node"] = node_role_counts["network_hub_node"]

    return df


def nodes_attributes_to_csv(G, filename):
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

    return df


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


def hex_to_rgba(value, alpha=1.0):
    value = value.lstrip('#')
    lv = len(value)
    ret = [int(value[i:i + lv // 3], 16)/256 for i in range(0, lv, lv // 3)]
    ret.append(alpha)
    return ret


def remap(x, oMin, oMax, nMin, nMax):
    #range check
    if oMin == oMax:
        print("Warning: Zero input range")
        return None

    if nMin == nMax:
        print("Warning: Zero output range")
        return None

    #check reversed input range
    reverseInput = False
    oldMin = min( oMin, oMax )
    oldMax = max( oMin, oMax )
    if not oldMin == oMin:
        reverseInput = True

    #check reversed output range
    reverseOutput = False   
    newMin = min( nMin, nMax )
    newMax = max( nMin, nMax )
    if not newMin == nMin :
        reverseOutput = True

    portion = (x-oldMin)*(newMax-newMin)/(oldMax-oldMin)
    if reverseInput:
        portion = (oldMax-x)*(newMax-newMin)/(oldMax-oldMin)

    result = portion + newMin
    if reverseOutput:
        result = newMax - portion
    
    return result


color_list = [
    "#d62728",    # red
    "#ffb0b0",    # light red
    "#ff00e6",    # pink
    "#ffb0f3",    # light pink
    "#ff7f0e",    # orange
    "#ffbb78",    # light orange
    "#8c564b",    # brown
    "#c49c94",    # light brown
    "#9467bd",    # purple
    "#c5b0d5",    # light purple
    "#bcbd22",    # mossy green 
    "#dbdb8d",    # light mossy green
    "#2ca02c",    # green
    "#98df8a",    # light green
    "#1f77b4",    # dark blue
    "#aec7e8",    # light blue
    "#17becf",    # sky blue
    "#9edae5",    # light sky blue
    "#7f7f7f",    # gray
    "#c7c7c7",    # light gray
]


color_rgba = [hex_to_rgba(c, alpha=1.0) for c in color_list]
cmap = ListedColormap(color_rgba, name='cmap')

color_rgba_transparent = [hex_to_rgba(c, alpha=0.8) for c in color_list]
cmap_transparent = ListedColormap(color_rgba_transparent, name='cmap_transparent')


special_colors = {
    "opitutus spp.": "black",                    # ???
    "eubacterium sp.": 'lime',                   # cellulose degrader
    "magnetospirillum sp.": 'lime',              # cellulose degrader
    "pleomorphomonas oryzae": "blue",            # nitrogen fixer
    "rhodopseudomonas palustris": "blue",        # nitrogen fixer
}


def visualize_network(G, layout="spring"):
    # list of 20 colors
    colors = []
    for i in G.nodes():
        colors.append(color_list[G.nodes[i]["modularity"]])

    vmin = 0
    vmax = len(color_list)

    plt.figure()

    if layout == "spring":
        # equivalent to fruchterman_reingold_layout
        nx.draw(G, pos=nx.spring_layout(G), node_color=colors)
    elif layout == "circular":
        nx.draw_circular(G, node_color=colors)
    else:
        nx.draw(G, node_color=colors)

    # add color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cb = plt.colorbar(sm)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.ax.invert_yaxis()

    plt.title(f"{G.name} {layout} layout")


def visualize_network_fancy(G):
    """
    Draw special nodes and regular nodes separately
    """
    vmin = 0
    vmax = len(color_list)
    
    # equivalent to fruchterman_reingold_layout
    layout = nx.spring_layout(G)
    
    special_node_list = []
    special_modularity = []
    special_eigenvector_centrality = []
    special_edge_color_list = []
    
    regular_node_list = []
    regular_modularity = []
    regular_eigenvector_centrality = []
    
    # list of 20 colors
    for i in G.nodes():
        if i in special_colors.keys():
            special_node_list.append(i)
            print(f"special {i} in module {G.nodes[i]['modularity']}")
            special_modularity.append(cmap(G.nodes[i]["modularity"]))
            special_eigenvector_centrality.append(G.nodes[i]["eigenvector_centrality"])
            special_edge_color_list.append(special_colors[i])
        else:
            regular_node_list.append(i)
            regular_modularity.append(cmap_transparent(G.nodes[i]["modularity"]))
            regular_eigenvector_centrality.append(G.nodes[i]["eigenvector_centrality"])
    
    # print("special_modularity:", special_modularity)
    
    # scale node sizes
    sp_sz_min = min(special_eigenvector_centrality)
    sp_sz_max = max(special_eigenvector_centrality)
    special_node_size = [remap(sz, sp_sz_min, sp_sz_max, 160, 600) for sz in special_eigenvector_centrality]
        
    re_sz_min = min(regular_eigenvector_centrality)
    re_sz_max = max(regular_eigenvector_centrality)
    regular_node_size = [remap(sz, re_sz_min, re_sz_max, 160, 600) for sz in regular_eigenvector_centrality]
    
    print(f"{G.name}: specials {len(special_edge_color_list)} | eigenvector centrality min: {min(sp_sz_min, re_sz_min)} | max: {min(sp_sz_max, re_sz_max)})")

    fig = plt.figure(figsize=(15, 15))

    nx.draw_networkx(G, with_labels=False, node_size=regular_node_size, \
        pos=layout,nodelist=regular_node_list, node_color=regular_modularity)
    special_nodes = nx.draw_networkx_nodes(G, node_size=special_node_size, \
        pos=layout,nodelist=special_node_list, node_color=special_modularity)

    # make special nodes have special outlines
    special_nodes.set_edgecolor(special_edge_color_list)
    
    # add color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cb = plt.colorbar(sm)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.ax.invert_yaxis()

    plt.title(f"{G.name} Network")


def get_leading_eigengenes(G, df_counts_rel_stan, set_node_module_membership=False):
    """ 
    G: networkx graph
    df_counts_rel_stan: pandas dataframe of relative abundance counts 
        standardized to have 0 mean and 1 variance
    """

    sample_names = list(df_counts_rel_stan.index)
    num_samples = len(sample_names)

    communities = nx_comm.greedy_modularity_communities(G)
    num_modules = len(communities)

    leading_eigenvalues = {}
    leading_eigenvalues_explained_variance = {}
    leading_eigenvector = {}

    df = pd.DataFrame()
    eigengenes = np.zeros((num_modules, num_samples), dtype=np.float64)

    node_module_membership_dict = {}
    for i, module_OTUs in enumerate(communities):
        df_X_b = df_counts_rel_stan.loc[:, list(module_OTUs)].T
        X_b = df_X_b.to_numpy(dtype=np.float64)
        u, s, vh = np.linalg.svd(X_b)

        # arrange the singular values in decreasing order
        idx = s.argsort()[::-1]
        eigenValues = s[idx]
        eigenVectors = vh[:,idx]

        # leading eigengene info
        df.loc[i, "leading_eigenvalue"] = eigenValues[0]
        df.loc[i, "leading_eigenvalues_percentage_explained_variance"] = np.square(eigenValues[0]) / np.sum(np.square(eigenValues)) * 100
        leading_eigengene = eigenVectors[:, 0].T
        eigengenes[i, :] = leading_eigengene

        # calculate module membership
        if (set_node_module_membership):            
            for OTU in module_OTUs:
                assert G.nodes[OTU]["modularity"] == i, f"{OTU} modularity {G.nodes[OTU]['modularity']} does not match {i}!"

                OTU_profile = df_counts_rel_stan.loc[:, OTU].values
                node_module_membership_dict[OTU] = np.corrcoef(OTU_profile, leading_eigengene)[0, 1] # pearson correlation
    
    if (set_node_module_membership):
        print(f"setting node_module_membership for {G.name} graph")
        nx.set_node_attributes(G, node_module_membership_dict, 'node_module_membership')

    df_eigengenes = pd.DataFrame(eigengenes, columns=sample_names)

    return df, df_eigengenes



    
















    


    