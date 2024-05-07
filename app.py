import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.neighbors import radius_neighbors_graph
import numpy as np
import standardise as se
from scipy.sparse.csgraph import connected_components
from collections import Counter
from scipy.spatial.distance import euclidean
from scipy.cluster.hierarchy import dendrogram, linkage

def create_graph(dataframe, radius, celltype=False):
    if not celltype:
        df = dataframe
    else:        
        df = dataframe[dataframe["celltype"] == celltype]
    coordinates_list = df[['nucleus.x', 'nucleus.y']].values
    info_list = df[['cell.ID', 'celltype']]
    graph = radius_neighbors_graph(coordinates_list, radius, mode='distance', metric='minkowski', p=2, metric_params=None, include_self=False)
    return graph, info_list

def main():
    st.title('Spatial Data')
    filename = st.sidebar.selectbox('Select a file', os.listdir("./if_data"))
    df = pd.read_csv(os.path.join("./if_data", filename))

    df_mapping = pd.read_csv("./IF1_phen_to_cell_mapping.csv")
    df['phenotype'] = df['phenotype'].apply(se.standardize_phenotype)
    merged_df = pd.merge(df, df_mapping, on='phenotype', how='left')

    # Plot the spatial data
    fig, ax = plt.subplots()
    sns.scatterplot(x="nucleus.x", y="nucleus.y", data=merged_df, hue="celltype")
    plt.title('Cell data')
    plt.xlabel("x_column")
    plt.ylabel("y_column")
    plt.legend(title="Cell type", loc='upper right')
    plt.grid(True)
    st.pyplot(fig)

    radius = st.sidebar.slider('Select radius:', min_value=1, max_value=50, value=50)
    min_bcell = st.sidebar.slider('Select the minimum amount of bcells in a component:', min_value=1, max_value=30, value=15)

    celltype = 'Bcell'
    # Create graphs
    graph_bcell, info_list_bcell = create_graph(merged_df, radius, celltype)
    graph, info_list = create_graph(merged_df, 20)

    
    # Find components, filter for components with > 15 bcells
    n_compontents, labels = connected_components(csgraph=graph_bcell, directed=False, return_labels=True)
    counts = np.bincount(labels)
    numbers_above_15 = counts[counts >= min_bcell]
    indexes_above_threshold = np.where(counts >= min_bcell)[0]
    if len(indexes_above_threshold) == 0:
        st.write("No components with 15 or more Bcells found.")
        for i in range(min_bcell):
            indexes_above_threshold = np.where(counts >= (min_bcell - i))[0]
            if len(indexes_above_threshold) != 0:
                st.write(f"{len(indexes_above_threshold)} component(s) with {min_bcell - i} Bcells found.")
                break
    else:
        st.write(f'{len(indexes_above_threshold)} component(s) found.')

    # Create a dictionary for each component's bcell amount
    bcell_component_count = {}
    for comp in indexes_above_threshold:
        bcell_component_count[comp] = counts[comp]

    # Find cell.ID and celltype for each bcell that's in a component with > 15 bcells
    components_cell_dict = {}
    for index in indexes_above_threshold:
        cell_dict = {}
        for cell in np.where(labels == index)[0]:
            cell_dict[int(cell)] = info_list_bcell.iloc[cell]
        components_cell_dict[int(index)] = cell_dict
    

    # Find the surroundings for each component
    component_neighbor_cells = {}
    for component_num in components_cell_dict.keys():
        neighbors = []
        for cell_num in components_cell_dict[component_num]:
            cell_id = components_cell_dict[component_num][cell_num]['cell.ID']
            row_number = info_list.index[info_list['cell.ID'] == cell_id][0]
            neighbors += [tuple(info_list.iloc[index]) for index in np.where(graph[row_number].toarray()[0] != 0)[0]]
        component_neighbor_cells[component_num] = list(set(neighbors))


    # Calculate the percentage of celltypes in each components surroundings
    component_celltype_percentages = {}
    for component_number in component_neighbor_cells.keys():
        data = component_neighbor_cells[component_number]
        celltypes = [celltype for _, celltype in data]
        counts = Counter(celltypes)
        total_cells = len(celltypes)
        percentages = {celltype: count / total_cells * 100 for celltype, count in counts.items()}
        component_celltype_percentages[component_number] = percentages
    
    components = list(component_celltype_percentages.keys())
    cell_types = set()
    if len(indexes_above_threshold) > 1:   
        # Pad with zeros
        for key in component_celltype_percentages:
            cell_types.update(component_celltype_percentages[key].keys())
        normalized_data = np.zeros((len(components), len(cell_types)))
        
        for i, key in enumerate(components):
            total = sum(component_celltype_percentages[key].values())
            for j, cell_type in enumerate(cell_types):
                normalized_data[i, j] = component_celltype_percentages[key].get(cell_type, 0) / total
        
        # Calculate the distance matrix for components
        distance_matrix = np.zeros((len(components), len(components)))
        for i in range(len(components)):
            for j in range(i, len(components)):
                distances = [euclidean(normalized_data[i], normalized_data[j])]
                distance_matrix[i, j] = distances[0]
                distance_matrix[j, i] = distances[0]
        
        Z = linkage(distance_matrix, 'ward')
        
        # Plot clustering dendrogram
        fig2 = plt.figure(figsize=(10, 5))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('Component')
        plt.ylabel('Distance')
        dendrogram(
            Z,
            labels=components,
            leaf_rotation=90.,
            leaf_font_size=12.,
        )
        st.pyplot(fig2)   
    

    # Plot histograms
    component_for_plotting = st.sidebar.selectbox('Select a component for plotting', components)
    component_celltype_percentages[component_for_plotting] = dict(sorted(component_celltype_percentages[component_for_plotting].items()))

    fig3 = plt.figure(figsize=(10, 5))
    plt.bar(component_celltype_percentages[component_for_plotting].keys(), component_celltype_percentages[component_for_plotting].values())
    plt.xlabel('Cell Type')
    plt.ylabel('Percentage of Cells')
    plt.title(f'Cell type percentages for component {component_for_plotting}')
    st.write(f'Bcell amount for component {component_for_plotting}: {bcell_component_count[component_for_plotting]}')
    st.pyplot(fig3)



if __name__ == "__main__":
    main()