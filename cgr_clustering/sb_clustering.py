from CGRtools.periodictable import Element, Core, At, DynamicElement
from CGRtools.containers import *
from CGRtools.containers.bonds import DynamicBond
from typing import Union, Dict, Tuple, List, Set, Optional
from IPython.display import display, SVG, HTML

from collections import defaultdict
import copy

from tqdm import tqdm

from CGRtools.algorithms.depict import *
import collections
from collections import defaultdict
from uuid import uuid4
import CGRtools # Assuming CGRTools is importable
import math # Import math for hypot if needed in original depict code

from functools import partial
from math import atan2, sin, cos, hypot

from typing import Dict, Any, Tuple

def extract_strat_bonds(target_cgr):
    """Extracts strategic bonds (order=None, p_order!=None)."""
    result = []
    seen = set()
    for atom1, bond_set in target_cgr._bonds.items():
        for atom2, bond in bond_set.items():
            if atom1 >= atom2:
                continue
            if bond.order is None and bond.p_order is not None:
                bond_key = tuple(sorted((atom1, atom2)))
                if bond_key not in seen:
                     seen.add(bond_key)
                     result.append(bond_key)
    return sorted(result)

def group_by_strat_bonds(data_dict, use_strat = True):
    """
    Groups rg_cgr objects based on their strategic bonds or 

    Args:
        data_dict: Dictionary mapping node_id to rg_cgr objects.

    Returns:
        Dictionary with groups keyed by '{length}.{index}' containing
        'sb_cgr', 'node_ids', and 'strat_bonds'.
    """
    temp_groups = collections.defaultdict(lambda: {'node_ids': [], 'sb_cgr': None, 'strat_bonds': None})

    # 1. Initial grouping based on the content of strategic bonds
    for node_id, rg_cgr in data_dict.items():
        strat_bonds_list = extract_strat_bonds(rg_cgr)
        if use_strat == True:
            group_key = tuple(strat_bonds_list)
        else:
            group_key = str(rg_cgr)

        if not temp_groups[group_key]['node_ids']: # First time seeing this group
            temp_groups[group_key]['sb_cgr'] = rg_cgr # Store the first CGR as representative
            temp_groups[group_key]['strat_bonds'] = strat_bonds_list # Store the actual list

        temp_groups[group_key]['node_ids'].append(node_id)
        temp_groups[group_key]['node_ids'].sort() # Keep node_ids sorted for consistency
        # temp_groups[group_key]['group_size'] = len(temp_groups[group_key]['node_ids'])
    for group_key in temp_groups.keys():
        temp_groups[group_key]['group_size'] = len(temp_groups[group_key]['node_ids'])

    # 2. Format the output dictionary with desired keys '{length}.{index}'
    final_grouped_results = {}
    group_indices = collections.defaultdict(int) # To track index for each length

    # Sort items by length of bonds first, then potentially by bonds themselves for consistent indexing
    # Sorting by the group_key (tuple of tuples) provides a deterministic order
    sorted_groups = sorted(temp_groups.items(), key=lambda item: (len(item[0]), item[0]))

    for group_key, group_data in sorted_groups:
        num_bonds = len(group_data['strat_bonds'])
        group_indices[num_bonds] += 1 # Increment index for this length (1-based)
        final_key = f"{num_bonds}.{group_indices[num_bonds]}"
        final_grouped_results[final_key] = group_data

    return final_grouped_results

def merge_groups(data, key1, key2):
    """
    Merge the group at key2 into the group at key1 within the provided dictionary.
    
    Parameters:
    - data (dict): Original dictionary containing the groups.
    - key1 (str): The key of the target group to preserve.
    - key2 (str): The key of the source group whose data will be merged, then deleted.
    
    Returns:
    - dict: A new dictionary with the merged result.
    """
    # Work on a shallow copy of the main dict to avoid mutating original
    merged = copy.deepcopy(data)
    
    # If key2 is not present or is empty, return the copy unchanged
    if not key2 or key2 not in merged:
        return merged
    
    group1 = merged[key1]
    group2 = merged[key2]
    
    # Merge node_ids
    group1['node_ids'].extend(group2['node_ids'])
    
    # Update group_size
    group1['group_size'] += group2['group_size']
    
    # Remove the old group
    del merged[key2]

    return merged

def fix_dict_key_order(d: Dict[str, Any]) -> Dict[str, Any]:
    """
    Renumber keys "M.N" so that for each M the Ns run 1,2,3… in the order seen.
    Prints any old→new key mappings.
    Returns a new dict with values preserved.
    """
    next_minor: Dict[str, int] = {}
    new_dict: Dict[str, Any] = {}
    changes: list[Tuple[str, str]] = []

    for old_key, value in d.items():
        # split major and ignore original minor
        major, _ = old_key.split('.', 1)
        # init counter if first time seeing this major
        if major not in next_minor:
            next_minor[major] = 1

        new_key = f"{major}.{next_minor[major]}"
        next_minor[major] += 1

        if new_key != old_key:
            changes.append((old_key, new_key))

        new_dict[new_key] = value

    # report what changed
    if changes:
        print("Renamed keys:")
        for old, new in changes:
            print(f"  {old} → {new}")

    return new_dict