import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd
import numpy as np

def get_gencode_maps(gencode_fn = "../data/external/gen_v26_mapping.csv"):
    """Get mapping from ENSG to gene symbol and vice versa
    """
    v26_mapping = pd.read_csv(gencode_fn)
    ens2sym = {v26_mapping['gene_id'].values[i].split('.')[0]:v26_mapping['gene_name'].values[i] for i in range(len(v26_mapping))}
    sym2ens = {j:i for i,j in ens2sym.items()}
    return(ens2sym, sym2ens)

def read_pandas(panda_fn):
    """Read pandas dataframes"""
    panda_data = {}
    for key, datafn in panda_fn.items():
        print(key)
        panda_data[key] = pd.read_csv(datafn, sep = ' ', index_col = 0)
    return(panda_data)

def read_dragons(dragon_fn):
    """Read dragon dataframes"""
    dragon_data = {}
    for key, datafn in dragon_fn.items():
        print(key)
        dragon_data[key] = pd.read_csv(datafn, sep = '\t', index_col = 0)
    return(dragon_data)

def get_indegree(panda_data, th = 0.5):
    """Get indegree (gene degree) for PANDA dataframes"""
    df = pd.DataFrame()
    for key, datafn in panda_data.items():
        print(key)
        temp = pd.DataFrame((datafn*1*(datafn>th)).sum(axis = 0))
        temp = temp.rename(columns={0:key})
        df = df.merge(temp, left_index=True, right_index=True, how = 'outer').fillna(-1)
    return(df)
        
def get_outdegree(panda_data, th = 0.5):
    """Get outdegree (TF degree) for PANDA dataframes"""
    df = pd.DataFrame()
    for key, datafn in panda_data.items():
        print(key)
        temp = pd.DataFrame((datafn*1*(datafn>th)).sum(axis = 1))
        temp = temp.rename(columns={0:key})
        df = df.merge(temp, left_index=True, right_index=True, how = 'outer').fillna(-1)
    return(df)
        
        
def get_long_dragon(dragons):
    """melt all dragons
    returns dataframe with melt dragons. 
    Pass a dictionary with dragons"""
    
    for k,key in enumerate(dragons.keys()):
        
        temp = dragons[key]
        # Step 1: Create a mask for the upper triangular matrix excluding the diagonal
        mask = np.triu(np.ones(temp.shape), k=1).astype(bool)

        # Step 2: Apply the mask to filter the DataFrame
        df_upper_tri = temp.where(mask)
        
        dragon = df_upper_tri.stack().reset_index()
        print(dragon)
        dragon.columns = ['node1', 'node2', key]
        print(dragon)
        if k==0:
            df = dragon
        else:
            df = df.merge(dragon, on=['node1', 'node2'], how = 'outer')
    return(df)

def draw_both_graphs(G, fig_name):
    
    """Draw both CMS2 and CMS4 graphs. saves the figure as fig_name
    """
    pos = nx.spring_layout(G, seed=200, iterations=35, k =  1/np.sqrt(len(G.nodes)))  # positions for all nodes
    #pos = nx.spectral_layout(G, scale=.1, center=None, dim=2)
    weights2 = nx.get_edge_attributes(G,'cms2').values()
    weights4 = nx.get_edge_attributes(G,'cms4').values()
    edge_labels = nx.get_edge_attributes(G,'edge_label_oi').values()
    edge_width = list(nx.get_edge_attributes(G,'edge_label_width').values())

    f,ax = plt.subplots(1,2,figsize = (15,10))
    options_node = {"node_size": 100, "alpha": 0.5}

    node_meth = [i  for i in list(G.nodes) if i.endswith('methylation') ]
    node_expr = [i  for i in list(G.nodes) if i.endswith('expression') ]

    # Methylation and expression nodes for CMS2
    nx.draw_networkx_nodes(G, pos, nodelist=node_meth, node_color="#66c2a5", **options_node, ax = ax[0])#, edgecolors='k')
    nx.draw_networkx_nodes(G, pos, nodelist=node_expr, node_color="#fc8d62", **options_node, ax = ax[0])#, edgecolors='k')
    
    # Methylation and expression nodes for CMS4
    nx.draw_networkx_nodes(G, pos, nodelist=node_meth, node_color="#66c2a5", **options_node, ax = ax[1])#, edgecolors='k')
    nx.draw_networkx_nodes(G, pos, nodelist=node_expr, node_color="#fc8d62", **options_node, ax = ax[1])#, edgecolors='k')

    # Node labels for methylation
    nx.draw_networkx_labels(G, pos, {node_meth[k]:node_meth[k].split('_')[0] for k in range(len(node_meth))}, font_size=8, ax = ax[0], verticalalignment='bottom', horizontalalignment='left')
    nx.draw_networkx_labels(G, pos, {node_meth[k]:node_meth[k].split('_')[0] for k in range(len(node_meth))}, font_size=8, ax = ax[1], verticalalignment='bottom', horizontalalignment='left')

    # Node labels for expression
    nx.draw_networkx_labels(G, pos, {node_expr[k]:node_expr[k].split('_')[0] for k in range(len(node_expr))}, font_size=8, ax = ax[0], verticalalignment='bottom', horizontalalignment='left')
    nx.draw_networkx_labels(G, pos, {node_expr[k]:node_expr[k].split('_')[0] for k in range(len(node_expr))}, font_size=8, ax = ax[1], verticalalignment='bottom', horizontalalignment='left')

    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.8, edge_color=[50*i for i in list(weights2)], edge_cmap=plt.cm.RdBu_r,  ax = ax[0], width = edge_width)
    nx.draw_networkx_edges(G, pos, alpha=0.8, edge_color=[50*i for i in list(weights4)], edge_cmap=plt.cm.RdBu_r, ax = ax[1], width = edge_width)

    #nx.draw_networkx_edge_labels(G, pos, edge_labels = nx.get_edge_attributes(G,'edge_label_oi'), ax = ax[0],  bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
    #nx.draw_networkx_edge_labels(G, pos, edge_labels = nx.get_edge_attributes(G,'edge_label_oi'), ax = ax[1],  bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})

    #nx.draw_networkx_labels(G_cms2, pos, {list(G_cms2.nodes)[k]:list(G_cms2.nodes)[k].split('_')[0] for k in range(len(G_cms2.nodes))}, font_size=9)


    meth_node = Line2D([0], [0], marker='o', color='#66c2a5', label='methylation',
                            markersize=10)
    expr_node = Line2D([0], [0], marker='o', color='#fc8d62', label='expression',
                            markersize=10)

    ax[0].legend(handles=[meth_node, expr_node])
    ax[1].legend(handles=[meth_node, expr_node])
    
    
    ax[0].set_title('CMS2')
    ax[1].set_title('CMS4')
    plt.tight_layout()

    f.savefig(fig_name, bbox_inches = 'tight')
    
    return(f)

