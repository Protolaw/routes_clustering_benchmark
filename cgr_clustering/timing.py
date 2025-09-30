import time
from tqdm import tqdm
from cgr_clustering.cgr import process_single_route, compose_all_sb_cgrs
from cgr_clustering.sb_clustering import group_by_strat_bonds
from cgr_clustering.aizynth_converter import (
    extract_pathway_aizynthfinder,
    route_smi_2_cgr
)
import numpy as np
from itertools import chain
from aizynthfinder.analysis.routes import RouteCollection

import matplotlib.pyplot as plt

def extract_route_cgr_sized(filtered_route_collection, size):
    route_cgrs_dict = {}
    for i, data in tqdm(enumerate(filtered_route_collection)):
        try:
            root = data.dicts[0]
        except:
            root = data['dict']
        pathway = extract_pathway_aizynthfinder(root)
        cgr_pathway = route_smi_2_cgr(pathway)
        route_cgr = process_single_route(cgr_pathway) 
        route_cgrs_dict[i] = route_cgr
        if len(route_cgrs_dict) > size:
            break
    return route_cgrs_dict

def sb_clustering_timed_run(filtered_route_collection, size, ):
    # Section 1
    t0 = time.time()
    route_cgrs_dict = extract_route_cgr_sized(filtered_route_collection, size)
    t1 = time.time()

    # Section 2
    sb_cgrs_dict = compose_all_sb_cgrs(route_cgrs_dict)
    t2 = time.time()

    # Section 3
    sbp_groups = group_by_strat_bonds(sb_cgrs_dict, use_strat=False)
    # cluster_arr = coll.cluster(len(sbp_groups), distances_model="ted")
    t3 = time.time()

    return (t1 - t0,      # extract_route_cgr
            t2 - t1,      # compose_all_sb_cgrs
            t3 - t2,      # clustering
            t3 - t0,      # total start→end
            )

def multitime_run(filtered_route_collection, sizes, ntries=5, method='sb'):
    times = {}
    stds = []
    for size in sizes:
        runs = [sb_clustering_timed_run(filtered_route_collection, size) for _ in range(ntries)]
        sec1, sec2, sec3, total = zip(*runs)
        arr = np.array(total)

        avg_sec1  = sum(sec1) / len(sec1)
        avg_sec2  = sum(sec2) / len(sec2)
        avg_sec3  = sum(sec3) / len(sec3)
        avg_total = sum(total) / len(total)

        print(f"Size {size}:")
        print(f"  • extract_route_cgr:    {avg_sec1:.4f} s")
        print(f"  • compose_all_sb_cgrs:   {avg_sec2:.4f} s")
        print(f"  • clustering:   {avg_sec3:.4f} s")
        print(f"  • total (start→end):     {avg_total:.4f} s")
        print("-" * 40)
        stds.append(arr.std(ddof=1))
        times[size] = {
            'sec1': avg_sec1,
            'sec2': avg_sec2,
            'sec3': avg_sec3,
            'total': avg_total
        }
    return times

def multitime_ted_cluster(filtered_route_collection, sizes, ntries=5):
    """
    Measure the time taken to cluster route collections of different sizes.
    
    Args:
        filtered_route_collection: The collection of routes to cluster.
        sizes: List of sizes to test.
        
    Returns:
        means: List of mean times for each size.
        stds: List of standard deviations for each size.
    """
    # To hold the per‐size statistics
    means = []
    stds = []

    for size in tqdm(sizes):
        secs = []
        for _ in range(ntries):
            x_route_collection = route_collection_for_cluster_size(filtered_route_collection, size=size)
            start = time.time()
            x_route_collection.cluster(17, distances_model="ted")
            end = time.time()
            secs.append(end - start)
        
        arr = np.array(secs)
        means.append(arr.mean())
        stds.append(arr.std(ddof=1))  # use ddof=1 for sample std, or leave default for population std
    print("Means:", means)
    print("Stds:", stds)
    return means, stds

def route_collection_for_cluster_size(route_collections, size=0):
    comb_filt_reaction_trees = [RouteCollection([i]) for i in route_collections.reaction_trees]
    all_trees = list(chain.from_iterable(rc.reaction_trees for rc in comb_filt_reaction_trees))[:size]

    # 2) Likewise, grab all of the MctsNode objects
    all_nodes = list(chain.from_iterable(rc.nodes for rc in comb_filt_reaction_trees))[:size]

    # 3) Now build your new, big collection:
    x_route_collection = RouteCollection(
        all_trees,
        nodes=all_nodes,
    )
    return x_route_collection


def plot_results(sizes, times_sizes_ted, times_sizes_sb):
    # Plot
    plt.figure()
    plt.plot(sizes, list(times_sizes_ted), marker='o', label='TED')
    plt.plot(sizes, list(times_sizes_sb), marker='s', color='red', label='SB-CGR total')

    # Labels and title
    plt.xlabel('Sample Size')
    plt.ylabel('Execution Time (s)')
    plt.title('Clustering Time vs. Sample Size')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)

    # Show the plot
    plt.show()