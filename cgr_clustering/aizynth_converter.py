# aizynth_converter.py

from aizynthfinder.analysis import TreeAnalysis
from aizynthfinder.analysis.routes import RouteCollection
from aizynthfinder.reactiontree import ReactionTree
from typing import List, Any, Dict
from itertools import chain


def extract_routes_from_tree(app):
    """Extracts all solved routes from the AiZynthFinder search tree."""
    tree_analysis = TreeAnalysis(app.finder.tree)
    nodes = list(tree_analysis.search_tree.graph())

    solved = [node for node in nodes if not node.children and node.state.is_solved]

    route_collections = []
    for node in solved:
        rt = node.to_reaction_tree()
        rc = RouteCollection([rt], nodes=[node])
        route_collections.append(rc)
    return route_collections


def filter_unique_routes(route_collections: RouteCollection) -> RouteCollection:
    """Filters a RouteCollection to keep only unique retrosynthetic routes."""

    all_trees = list(chain.from_iterable(rc.reaction_trees for rc in route_collections))
    all_nodes = list(chain.from_iterable(rc.nodes for rc in route_collections))
    
    single_route_collection = RouteCollection(
        all_trees,
        nodes=all_nodes,
    )

    seen_hashes = set()
    unique_trees: List[ReactionTree] = []
    unique_nodes: List[Any] = []
    unique_dicts: List[Dict[str, Any]] = []

    for tree, node, dct in zip(
        single_route_collection.reaction_trees,
        single_route_collection.nodes,
        single_route_collection.dicts
    ):
        h = tree.hash_key()
        if h in seen_hashes:
            continue
        seen_hashes.add(h)
        unique_trees.append(tree)
        unique_nodes.append(node)
        unique_dicts.append(dct)

    return RouteCollection(
        unique_trees,
        nodes=unique_nodes,
        dicts=unique_dicts
    )


def extract_pathway_aizynthfinder(node, parent_smiles=None):
    """Recursively extracts a pathway from a reaction tree node."""
    pathway = []
    if node.get('type') == 'reaction':
        for child in node.get('children', []):
            if child.get('type') == 'mol' and 'children' in child:
                for sub in child['children']:
                    pathway.extend(extract_pathway_aizynthfinder(sub, child['smiles']))
        reactants = [c['smiles'] for c in node['children'] if c['type']=='mol'][::-1]
        pathway.append([reactants, parent_smiles])
    else:
        for child in node.get('children', []):
            pathway.extend(extract_pathway_aizynthfinder(child, node.get('smiles')))
    return pathway
