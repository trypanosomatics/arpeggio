"""
Tests for the mmCIF readers in arpeggio.core.protein_reader.

Fixtures (made-up coordinates): a sulfate numbered 401 in chain B, followed by a copy numbered 401 in
chain BA, as in a biological assembly built from symmetry copies. One file has alternate locations A/B
(occupancy 0.5 each), the other none.

Run from the repository root:  python -m unittest discover -s arpeggio/tests
"""

import os
import unittest

from openbabel import openbabel as ob
from Bio.PDB.Atom import DisorderedAtom

from arpeggio.core import protein_reader

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ALTLOC = os.path.join(DATA, "chain_boundary_altloc.cif")
NO_ALTLOC = os.path.join(DATA, "chain_boundary_no_altloc.cif")


def biopython_serials(structure):
    """
    Serial numbers of all atoms, alternate locations included (as InteractionComplex collects them).
    """
    serials = []
    for atom in structure.get_atoms():
        children = atom.child_dict.values() if isinstance(atom, DisorderedAtom) else [atom]
        serials.extend(child.serial_number for child in children)
    return sorted(serials)


def openbabel_serials(mol):
    return sorted(atom.GetId() for atom in ob.OBMolAtomIter(mol))


class TestChainBoundary(unittest.TestCase):

    def test_altloc_copies_keep_their_own_chains(self):
        structure = protein_reader.read_mmcif_to_biopython(ALTLOC)
        chains = {chain.id: chain for chain in structure[0]}
        self.assertEqual(sorted(chains), ["B", "BA"])
        for chain_id, first_serial in (("B", 1), ("BA", 11)):
            residue = chains[chain_id][("H_SO4", 401, " ")]
            serials = sorted(
                child.serial_number for atom in residue for child in atom.child_dict.values()
            )
            self.assertEqual(serials, list(range(first_serial, first_serial + 10)))

    def test_altloc_biopython_and_openbabel_serials_match(self):
        # a mismatch here is what makes InteractionComplex raise OBBioMatchError
        bio = biopython_serials(protein_reader.read_mmcif_to_biopython(ALTLOC))
        obabel = openbabel_serials(protein_reader.read_mmcif_to_openbabel(ALTLOC))
        self.assertEqual(bio, list(range(1, 21)))
        self.assertEqual(bio, obabel)

    def test_openbabel_residues_follow_chains(self):
        mol = protein_reader.read_mmcif_to_openbabel(ALTLOC)
        residues = [(res.GetChainNum(), res.GetNumAtoms()) for res in ob.OBResidueIter(mol)]
        self.assertEqual(residues, [(1, 10), (2, 10)])

    def test_without_altlocs(self):
        structure = protein_reader.read_mmcif_to_biopython(NO_ALTLOC)
        chains = {chain.id: chain for chain in structure[0]}
        self.assertEqual(sorted(chains), ["B", "BA"])
        for chain_id in ("B", "BA"):
            self.assertEqual(len(chains[chain_id][("H_SO4", 401, " ")]), 5)
        bio = biopython_serials(structure)
        self.assertEqual(bio, openbabel_serials(protein_reader.read_mmcif_to_openbabel(NO_ALTLOC)))


if __name__ == "__main__":
    unittest.main()
