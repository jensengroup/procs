import sys
import numpy
from Bio.PDB import *

def get_chi(residue):
        if residue.get_resname() == 'ALA':
                return []
        if residue.get_resname() == 'GLY':
                return []
        if residue.get_resname() == 'ARG':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['NE'].get_vector()
                sc_atom7 = residue['CZ'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom6, sc_atom7)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'ASN':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['OD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'ASP':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['OD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'CYS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['SG'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'GLU':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['OE1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'GLN':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['OE1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'HIS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD2'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'ILE':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG1'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'LEU':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'LYS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['CE'].get_vector()
                sc_atom7 = residue['NZ'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom6, sc_atom7)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'MET':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['SD'].get_vector()
                sc_atom6 = residue['CE'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'PHE':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'PRO':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom1)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom1, sc_atom2)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'SER':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['OG'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'THR':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['OG1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'VAL':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        else:
                return "Error reading side chain angles"

def get_phi_psi_angles(chain):

    polypeptide = Polypeptide.Polypeptide(chain)

    phi_psi_disordered_list = polypeptide.get_phi_psi_list()

    phi_psi_ordered_dictionary = dict()

    for res_index, residue in enumerate(polypeptide):
        phi_psi_ordered_dictionary[residue.id[1]] = phi_psi_disordered_list[res_index]

    print phi_psi_ordered_dictionary

    return phi_psi_ordered_dictionary
    



def load_chain(filename):
    for model in PDBParser().get_structure("chain", "5PTI.pdb") :
        for chain in model:
            return chain



filename = sys.argv[1]

chain = load_chain(filename)

# Calculate phi/psi angles once
print get_phi_psi_angles(chain)

for residue in chain:
    print residue



