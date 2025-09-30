import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter

import seaborn as sns

def plot_distribution(group_sizes, clusters=None, method='LSTM'):
    if len(group_sizes) > 10:

        # Count how often each value appears
        freq = Counter(group_sizes)
        values = sorted(freq.keys())
        counts = [freq[v] for v in values]

        # Plot
        plt.figure(figsize=(9, 6))
        plt.bar(values, counts, width=0.8)
        plt.xlabel('Size of Clusters (Number of routes)')
        plt.ylabel('Frequency')
        plt.title(f'Distribution of {len(group_sizes)} Clusters')
        plt.xticks(values)
        # plt.grid(axis='y', linestyle='--', linewidth=0.5)
        plt.grid(False)

        # Add count labels above each bar
        for x, y in zip(values, counts):
            plt.text(x, y + max(counts)*0.01, str(y), ha='center', va='bottom')

        plt.tight_layout()

    else:
        sns.set_style("whitegrid")
        labels = clusters.keys()
        total_num = sum(group_sizes)
        fig, ax = plt.subplots(figsize=(8, 8))
        wedges, texts, autotexts = ax.pie(
            group_sizes,
            autopct='%1.1f%%',
            startangle=140,
            pctdistance=1.1,
            colors=sns.color_palette("pastel")
        )
        plt.setp(autotexts, size=13, weight="bold")
        ax.legend(wedges, labels, title=f"Clusters", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=13, title_fontsize=13)
        ax.set_title(f'Cluster size distribution - {method} for {total_num} routes')

    plt.show()