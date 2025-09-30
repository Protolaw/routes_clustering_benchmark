from chython import ReactionContainer as ReactionContainerChython
from chython import smiles as smiles_chython

from CGRtools import smiles as smiles_cgrtools
from CGRtools.containers import ReactionContainer, CGRContainer
from CGRtools.containers.bonds import DynamicBond
from CGRtools import smiles as smiles_cgrtools
from CGRtools.algorithms.depict import *
from collections import defaultdict

from functools import partial
from math import atan2, sin, cos, hypot


def route_smi_2_cgr(pathway, reverse=False): # True for AiZynthFInder, False for ASKCOS
    """Converts a pathway of SMILES strings to a list of CGRs."""
    cgr_pathway = []
    inversed_pathway = pathway[::-1] if reverse else pathway
    for reaction_str in inversed_pathway:
        reactants = []
        product = smiles_chython(reaction_str[1])
        for reactant_smiles in reaction_str[0]:
            reactant = smiles_chython(reactant_smiles)
            reactant.kekule()
            reactant.implicify_hydrogens()
            reactant.thiele()
            reactants.append(reactant)
        reaction = ReactionContainerChython(reactants=reactants, products = [product])
        reaction.reset_mapping(keep_reactants_numbering=False)
        reaction_cgrtools = smiles_cgrtools(format(reaction, "m"))
        cgr_pathway.append(reaction_cgrtools)
    return cgr_pathway

def find_remap(lst):
    """
    Given a sorted list `lst` whose true length N is known to be len(lst),
    returns a dict mapping each value > N in lst to the missing values in 1..N.

    Example:
      L = [1,2,...,18,20,21,22,23]  # len=22
      => missing = [19]
         out_of_range = [23]
      => {23: 19}
    """
    N = len(lst)
    # 1) which values in the “ideal” 1..N are missing?
    missing = sorted(set(range(1, N+1)) - set(x for x in lst if x <= N))
    # 2) which values in lst have “overflowed” past N?
    out_of_range = sorted(x for x in lst if x > N)

    if len(missing) != len(out_of_range):
        raise ValueError(f"got {len(missing)} missing slots but {len(out_of_range)} overflow values")

    # 3) pair them up in ascending order
    return dict(zip(out_of_range, missing))

def process_single_route(cgr_pathway):
    for i, reaction in enumerate(cgr_pathway):
        if i == 0:
            cgr = reaction.compose()
            atoms = reaction.products[0]._atoms.keys()
            if reaction.products[0].atoms_count != max(atoms):
                remapper = find_remap(list(atoms))
                temp_num = max(cgr._atoms)+1
                for key, value in remapper.items():
                    save_val = int(value)
                    cgr.remap({value: temp_num, key: value, value: key})
        else:
            curr_product = reaction.products[0]
            curr_product.kekule()
            curr_product.implicify_hydrogens()
            curr_product.thiele()
            for reactant in decomposed.reactants:
                reactant.kekule()
                reactant.implicify_hydrogens()
                reactant.thiele()
                try:
                    if len(reactant) == len(curr_product):
                        curr_remap = next(curr_product.get_mapping(reactant))
                        curr_cgr = reaction.compose()
                        max_num = max(cgr._atoms) + 1
                        curr_decomposed = ReactionContainer.from_cgr(curr_cgr)
                        lg_remap = {}
                        for product in curr_decomposed.products:
                            curr_max_num = max(curr_cgr._atoms) + 1
                            if curr_max_num > max_num:
                                max_num = curr_max_num
                            if len(product) == len(curr_product):
                                continue
                            else:
                                for atom_num in product:
                                    lg_remap[atom_num] = max_num
                                    max_num += 1
                        curr_cgr.remap(lg_remap)
                        curr_cgr.remap(curr_remap)
                        new_reaction = ReactionContainer.from_cgr(curr_cgr)
                        cgr = curr_cgr.compose(cgr)
                except:
                    pass
        decomposed = ReactionContainer.from_cgr(cgr)
    target_cgr = [cgr.substructure(c) for c in cgr.connected_components][0]
    target_cgr = cgr_enhance(target_cgr)
    return target_cgr

def compose_sb_cgr(route_cgr: CGRContainer):
    """
    Reduces a Routes Condensed Graph of reaction (RouteCGR) by performing the following steps:

    1. Extracts substructures corresponding to connected components from the input RouteCGR.
    2. Selects the first substructure as the target to work on.
    3. Iterates over all bonds in the target RouteCGR:
       - If a bond is identified as a "leaving group" (its primary order is None while its original order is defined),
         the bond is removed.
       - If a bond has a modified order (both primary and original orders are integers) and the primary order is less than the original,
         the bond is deleted and then re-added with a new dynamic bond using the primary order (this updates the bond to the reduced form).
    4. After bond modifications, re-extracts the substructure from the target RouteCGR (now called the reduced RouteCGR or ReducedRouteCGR).
    5. If the charge distributions (_p_charges vs. _charges) differ, neutralizes the charges by setting them to zero.

    Args:
        route_cgr: The input RouteCGR object to be reduced.

    Returns:
        The reduced RouteCGR object.
    """
    # Get all connected components of the RouteCGR as separate substructures.
    cgr_prods = [route_cgr.substructure(c) for c in route_cgr.connected_components]
    target_cgr = cgr_prods[0]

    # Iterate over each bond in the target RouteCGR.
    bond_items = list(target_cgr._bonds.items())
    for atom1, bond_set in bond_items:
        bond_set_items = list(bond_set.items())
        for atom2, bond in bond_set_items:

            # Removing bonds corresponding to leaving groups:
            # If product bond order is None (indicating a leaving group) but an original bond order exists,
            # delete the bond.
            if bond.p_order is None and bond.order is not None:
                target_cgr.delete_bond(atom1, atom2)

            # For bonds that have been modified (not leaving groups) where the new (primary) order is less than the original:
            # Remove the bond and re-add it using the DynamicBond with the primary order for both bond orders.
            elif (
                type(bond.p_order) is int
                and type(bond.order) is int
                and bond.p_order != bond.order
            ):
                p_order = int(bond.p_order)
                target_cgr.delete_bond(atom1, atom2)
                target_cgr.add_bond(atom1, atom2, DynamicBond(p_order, p_order))

    # After modifying bonds, extract the reduced RouteCGR from the target's connected components.
    sb_cgr = [target_cgr.substructure(c) for c in target_cgr.connected_components][0]

    # Neutralize charges if the primary charges and current charges differ.
    if sb_cgr._p_charges != sb_cgr._charges:
        for num, charge in sb_cgr._charges.items():
            if charge != 0:
                sb_cgr._atoms[num].charge = 0
    sb_cgr = cgr_enhance(sb_cgr)
    return sb_cgr


def compose_all_sb_cgrs(route_cgrs_dict: dict):
    """
    Processes a collection (dictionary) of RouteCGRs to generate their reduced forms (ReducedRouteCGRs).

    Iterates over each RouteCGR in the provided dictionary and applies the compose_sb_cgr function.

    Args:
        route_cgrs_dict (dict): A dictionary where keys are identifiers (e.g., route numbers)
                                and values are RouteCGR objects.

    Returns:
        dict: A dictionary where each key corresponds to the original identifier from
              `route_cgrs_dict` and the value is the corresponding ReducedRouteCGR object.
    """
    all_sb_cgrs = dict()
    for num, cgr in route_cgrs_dict.items():
        all_sb_cgrs[num] = compose_sb_cgr(cgr)
    return all_sb_cgrs


def rotate_vector(x1, y1, x2, y2):
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 - sin_rad * y1, sin_rad * x1 + cos_rad * y1

class WideFormedBondDepictCGR(DepictCGR):
    """
    Like DepictCGR, but all bonds with order==None
    are drawn 4× wider than the standard bond width.
    """
    __slots__ = ()

    def _render_bonds(self):
        plane = self._plane
        config = self._render_config

        # get the normal width (default 1.0) and compute a 4× wide stroke
        normal_width = config.get('bond_width', 0.02) 
        wide_width = normal_width * 2.5

        broken = config['broken_color']
        formed = config['formed_color']
        dash1, dash2 = config['dashes']
        double_space = config['double_space']
        triple_space = config['triple_space']

        svg = []
        ar_bond_colors = defaultdict(dict)

        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            # invert Y for SVG
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order is None:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 4:
                if p_order == 4:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 1:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 3:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order is None:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed} stroke-width="{wide_width:.2f}""/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" 'f'stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{broken}"/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    # dx = dx // 1.4
                    # dy = dy // 1.4
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
            else:
                if p_order == 8:
                    svg.append(f'        <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}" stroke="{formed}" stroke-width="{wide_width:.2f}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')

        # aromatic rings - unchanged
        for ring in self.aromatic_rings:
            cx = sum(plane[x][0] for x in ring) / len(ring)
            cy = sum(plane[x][1] for x in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny = plane[n]
                mx, my = plane[m]
                aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
                if aromatic:
                    svg.append(aromatic)

            n, m = ring[-1], ring[0]
            nx, ny = plane[n]
            mx, my = plane[m]
            aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y, color):
        config = self._render_config

        dash1, dash2 = config['dashes']
        dash3, dash4 = config['aromatic_dashes']
        aromatic_space = config['cgr_aromatic_space']

        normal_width = config.get('bond_width', 0.02) 
        wide_width = normal_width * 2
        
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, -mn_y)

        if cr_y and aromatic_space / cr_y < .65:
            if cr_y > 0:
                r_y = aromatic_space
            else:
                r_y = -aromatic_space
                cr_y = -cr_y

            ar_x = aromatic_space * cr_x / cr_y
            br_x = mr_x - aromatic_space * (mr_x - cr_x) / cr_y

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, r_y, mn_x, mn_y)
            bn_x, bn_y = rotate_vector(br_x, r_y, mn_x, mn_y)
            
            if color:
                # print('color')
                return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}" x2="{bn_x + n_x:.2f}" ' \
                       f'y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{color}" stroke-width="{wide_width:.2f}"/>'
            elif color is None:
                dash3, dash4 = dash1, dash2
            return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}"' \
                   f' x2="{bn_x + n_x:.2f}" y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>'

def cgr_enhance(cgr: CGRContainer) -> str:
    
    CGRContainer._CGRContainer__render_aromatic_bond = (
        WideFormedBondDepictCGR._WideFormedBondDepictCGR__render_aromatic_bond
    )
    CGRContainer._render_bonds = WideFormedBondDepictCGR._render_bonds
    # CGRContainer.__render_aromatic_bond = WideFormedBondDepictCGR.__render_aromatic_bond
    CGRContainer._WideFormedBondDepictCGR__render_aromatic_bond = (
    WideFormedBondDepictCGR._WideFormedBondDepictCGR__render_aromatic_bond
    )
    cgr.clean2d()

    return cgr