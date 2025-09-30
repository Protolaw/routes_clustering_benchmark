import collections

def extract_clusters(data):
    clusters = {}
    for i, route in enumerate(data):
        print(i, route['graph']['cluster_id'])
        cluster_id = route['graph']['cluster_id'] + 2 # Adjusting to 1-based index
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(i)

    clusters = collections.OrderedDict(sorted(clusters.items()))
    return clusters

def parse_reaction_smiles(rxn_smiles):
    """
    Split a reaction SMILES into (reactants, product), where reactants is
    always a list and product is a string if there’s only one.
    """
    reactants_str, products_str = rxn_smiles.split('>>')
    reactants = reactants_str.split('.') if reactants_str else []
    products  = products_str.split('.')  if products_str  else []
    # unwrap single product
    product = products[0] if len(products) == 1 else products
    return reactants, product

def extract_pathway(nodes, parent_smiles=None):
    pathway = []
    for node in nodes:
        if node['type'] == 'reaction':
            rxn_smiles = node['smiles']
            reactants, product = parse_reaction_smiles(rxn_smiles)
            pathway.append((reactants, product))
            # print(node['smiles'])
    return pathway
