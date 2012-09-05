import numpy
from Bio.PDB import *


side_chain_atoms = dict()

side_chain_atoms['ALA'] = []
side_chain_atoms['GLY'] = []
side_chain_atoms['ARG'] = ['N',  'CA',  'CB',  'CG',  'CD',  'NE',  'CZ']
side_chain_atoms['ASN'] = ['N',  'CA',  'CB',  'CG',  'OD1']
side_chain_atoms['ASP'] = ['N',  'CA',  'CB',  'CG',  'OD1']
side_chain_atoms['CYS'] = ['N',  'CA',  'CB',  'SG']
side_chain_atoms['GLU'] = ['N',  'CA',  'CB',  'CG',  'CD',  'OE1']
side_chain_atoms['GLN'] = ['N',  'CA',  'CB',  'CG',  'CD',  'OE1']
side_chain_atoms['HIS'] = ['N',  'CA',  'CB',  'CG',  'CD2']
side_chain_atoms['ILE'] = ['N',  'CA',  'CB',  'CG1', 'CD1']
side_chain_atoms['LEU'] = ['N',  'CA',  'CB',  'CG', 'CD1']
side_chain_atoms['LYS'] = ['N',  'CA',  'CB',  'CG',  'CD',  'CE',  'NZ']
side_chain_atoms['MET'] = ['N',  'CA',  'CB',  'CG',  'SD',  'CE']
side_chain_atoms['PHE'] = ['N',  'CA',  'CB',  'CG',  'CD1']
side_chain_atoms['PRO'] = ['N',  'CA',  'CB',  'CG',  'CD']
side_chain_atoms['SER'] = ['N',  'CA',  'CB',  'OG']
side_chain_atoms['THR'] = ['N',  'CA',  'CB',  'OG1']
side_chain_atoms['TYR'] = ['N',  'CA',  'CB',  'CG', 'CD1']
side_chain_atoms['TRP'] = ['N',  'CA',  'CB',  'CG', 'CD1']
side_chain_atoms['VAL'] = ['N',  'CA',  'CB',  'CG1']


# Return a chi angle for a single residue
def get_chi(residue):

    chi = []

    atoms = side_chain_atoms[residue.get_resname()]

    if len(atoms) < 4:
	return chi

    for i in range(len(atoms) - 3):
        try:
            sc_atom0 = residue[atoms[0 + i]].get_vector()
            sc_atom1 = residue[atoms[1 + i]].get_vector()
            sc_atom2 = residue[atoms[2 + i]].get_vector()
            sc_atom3 = residue[atoms[3 + i]].get_vector()
    
            this_chi = calc_dihedral(sc_atom0,
                                     sc_atom1,
                                     sc_atom2,
                                     sc_atom3)
        except:
            print "Error: Missing atoms in side chain of ", residue.get_resname(), residue.id[1]
            this_chi = None

        chi.append(this_chi)

    return chi


# Return a dictionary of lists of phi/psi-angles for each residues
def get_phi_psi_angles(chain):

    polypeptide = Polypeptide.Polypeptide(chain)
    phi_psi_disordered_list = polypeptide.get_phi_psi_list()
    phi_psi_dictionary = dict()

    for res_index, residue in enumerate(polypeptide):

        phi_psi_tuple = phi_psi_disordered_list[res_index]
        phi_psi_dictionary[residue.id[1]] = [phi_psi_tuple[0], phi_psi_tuple[1]]

    return phi_psi_dictionary
    

# Return a dictionary of lists of chi-angles for each residues
def get_chi_angles(chain):

    chi_dictionary = dict()

    for residue in chain:
        chi_dictionary[residue.id[1]] = get_chi(residue)

    return chi_dictionary
	

# Load first chain of first model found in the pdb-file
def load_chain(filename):
    for model in PDBParser().get_structure("chain", filename) :
        for chain in model:
            return chain




