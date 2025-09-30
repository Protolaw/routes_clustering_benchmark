
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mutual_info_score, adjusted_rand_score, normalized_mutual_info_score
from scipy.stats import entropy
from itertools import product
import matplotlib.patches as patches
from scipy.optimize import linear_sum_assignment

def _extract_labels(first_dict, second_dict):
    """
    Internal helper to extract label maps and common items from two cluster dicts.
    """
    labels1 = {node: cid for cid, data in first_dict.items() for node in data['node_ids']}
    labels2 = {}
    for cid, nodes in second_dict.items():
        # second_dict may map to a list or to a dict with 'node_ids'
        if isinstance(nodes, dict) and 'node_ids' in nodes:
            node_list = nodes['node_ids']
        else:
            node_list = nodes
        for node in node_list:
            labels2[node] = cid
    common = sorted(set(labels1) & set(labels2))
    return labels1, labels2, common


def calculate_similarity_metrics(first_dict, second_dict):
    """
    Compute ARI, NMI, entropy, MI, and VI between two clusterings.
    Returns y_true and y_pred lists aligned on the common item set.
    """
    labels1, labels2, common = _extract_labels(first_dict, second_dict)
    if not common:
        print("No common items found between the two clusterings.")
        return [], []

    y_true = [labels1[n] for n in common]
    y_pred = [labels2[n] for n in common]

    # External metrics
    ari = adjusted_rand_score(y_true, y_pred)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    print(f"Adjusted Rand Index (ARI): {ari:.4f}")
    print(f"Normalized Mutual Information (NMI): {nmi:.4f}")

    # Info-theoretic metrics
    counts1 = np.unique(y_true, return_counts=True)[1]
    counts2 = np.unique(y_pred, return_counts=True)[1]
    H1 = entropy(counts1)
    H2 = entropy(counts2)
    MI = mutual_info_score(y_true, y_pred)
    VI = H1 + H2 - 2 * MI
    print(f"Entropy first: {H1:.4f}")
    print(f"Entropy second: {H2:.4f}")
    print(f"Mutual Information: {MI:.4f}")
    print(f"Variation of Information (VI): {VI:.4f}\n")

    return y_true, y_pred


def match_clusters(first_dict, second_dict):
    """
    Match clusters by maximizing average Jaccard similarity using the Hungarian algorithm.
    Returns list of matched tuples (cid1, cid2, score) and the average score.
    """
    labels1, labels2, common = _extract_labels(first_dict, second_dict)

    # Build set-based clusters restricted to common items
    clusters1 = {cid: set(data['node_ids']) & set(common)
                 for cid, data in first_dict.items()}
    clusters2 = {cid: set(nodes) & set(common)
                 for cid, nodes in second_dict.items()}

    keys1, keys2 = list(clusters1), list(clusters2)
    J = np.zeros((len(keys1), len(keys2)))
    for i, c1 in enumerate(keys1):
        for j, c2 in enumerate(keys2):
            s1, s2 = clusters1[c1], clusters2[c2]
            if s1 or s2:
                J[i, j] = len(s1 & s2) / len(s1 | s2)

    # Hungarian: maximize J by minimizing (1-J)
    row_ind, col_ind = linear_sum_assignment(1 - J)
    matched = [(keys1[i], keys2[j], J[i, j]) for i, j in zip(row_ind, col_ind)]
    avg_jaccard = float(np.mean([score for *_, score in matched]))
    return matched, avg_jaccard


def compute_contingency_matrix(y_true, y_pred):
    """
    Create a pandas contingency table from two label lists.
    """
    return pd.crosstab(
        pd.Series(y_true, name='first_partition'),
        pd.Series(y_pred, name='second_partition')
    )


def contingency_nodes(first_dict, second_dict, row_key, col_key):
    """Return sorted list of nodes shared by two specific clusters."""
    s1 = set(first_dict[row_key]['node_ids'])
    s2 = set(second_dict[col_key])
    return sorted(s1 & s2)


def display_contingency_matrix(first_dict, second_dict,
                               matched, cm, method='TED'):
    """Plot a heatmap of the contingency matrix, highlighting matched pairs."""
    labels_a = sorted(first_dict)
    labels_b = sorted(second_dict)
    np_cm = cm.to_numpy()

    # Determine highlight positions
    highlight = [(labels_a.index(r), c, score)
                 for r, c, score in matched]
    print(highlight)

    fig, ax = plt.subplots(figsize=(6 if method=='TED' else 12, 5 if method=='TED' else 3))
    im = ax.imshow(cm, aspect='auto', cmap='Greys',
                   vmin=0, vmax=np_cm.max()+1)

    # Annotate cells
    for i in range(np_cm.shape[0]):
        for j in range(np_cm.shape[1]):
            if np_cm[i, j] > 0:
                ax.text(j, i, np_cm[i, j], ha='center', va='center',
                        color='grey')

    # Add red highlights for matched pairs
    for i, j, score in highlight:
        if np_cm[i, j-1] == 0:
            continue
        rect = patches.Rectangle((j - 1.5, i - 0.5), 1, 1,
                                linewidth=2, edgecolor='red', facecolor='none',
                                alpha=0.5 + 0.5 * score)
        ax.add_patch(rect)
        ax.text(j-1, i, np_cm[i, j-1], ha='center', va='center',
                color='red', fontweight='bold')

    # Formatting
    ax.set_title(f"Contingency Matrix ({method} vs SB-CGR)")
    ax.set_xlabel(f"{method} Clusters")
    ax.set_ylabel("SB-CGR Clusters")
    ax.set_xticks(range(len(labels_b)))
    ax.set_xticklabels(labels_b)
    ax.set_yticks(range(len(labels_a)))
    ax.set_yticklabels(labels_a)
    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label('Number of Routes', rotation=270, labelpad=15)
    plt.tight_layout()
    plt.show()