import numpy
from Bio.PDB import *


# Atoms in the chain that forms chi1, chi2, etc. dihedral angles
SIDE_CHAINS = {'ALA' : [],
               'GLY' : [],
               'ARG' : ['N',  'CA',  'CB',  'CG',  'CD',  'NE',  'CZ'],
               'ASN' : ['N',  'CA',  'CB',  'CG',  'OD1'],
               'ASP' : ['N',  'CA',  'CB',  'CG',  'OD1'],
               'CYS' : ['N',  'CA',  'CB',  'SG'],
               'GLU' : ['N',  'CA',  'CB',  'CG',  'CD',  'OE1'],
               'GLN' : ['N',  'CA',  'CB',  'CG',  'CD',  'OE1'],
               'HIS' : ['N',  'CA',  'CB',  'CG',  'CD2'],
               'ILE' : ['N',  'CA',  'CB',  'CG1', 'CD1'],
               'LEU' : ['N',  'CA',  'CB',  'CG', 'CD1'],
               'LYS' : ['N',  'CA',  'CB',  'CG',  'CD',  'CE',  'NZ'],
               'MET' : ['N',  'CA',  'CB',  'CG',  'SD',  'CE'],
               'PHE' : ['N',  'CA',  'CB',  'CG',  'CD1'],
               'PRO' : ['N',  'CA',  'CB',  'CG',  'CD'],
               'SER' : ['N',  'CA',  'CB',  'OG'],
               'THR' : ['N',  'CA',  'CB',  'OG1'],
               'TYR' : ['N',  'CA',  'CB',  'CG', 'CD1'],
               'TRP' : ['N',  'CA',  'CB',  'CG', 'CD1'],
               'VAL' : ['N',  'CA',  'CB',  'CG1']}

# Atoms that form aromatic rings.
SIDE_CHAIN_RINGS = {'HIS' : [['CG',  'ND1',  'CD2',  'CE1',  'NE2']],
                    'TYR' : [['CG',  'CD1',  'CD2',  'CE1',  'CE2',  'CZ']],
                    'PHE' : [['CG',  'CD1',  'CD2',  'CE1',  'CE2',  'CZ']],
                    'TRP' : [['CG',  'CD1',  'CD2',  'NE1',  'CE2'],
                             ['CD2', 'CE2',  'CE3',  'CZ2',  'CZ3',  'CH2']]}

# Return a chi angle for a single residue
def calculate_chi_angles(residue):

    atoms = SIDE_CHAINS[residue.get_resname()]

    # For Ala and Gly (no chi angles)
    if len(atoms) < 4:
	return []

    chi = []

    for i in range(len(atoms) - 3):
        try:
            this_chi = calc_dihedral(residue[atoms[0 + i]].get_vector(),
                                     residue[atoms[1 + i]].get_vector(),
                                     residue[atoms[2 + i]].get_vector(),
                                     residue[atoms[3 + i]].get_vector())
        except:
            print "Error: Missing atoms in side chain of ", 
            print residue.get_resname(), residue.id[1]

            this_chi = None

        chi.append(this_chi)

    return chi


# Return a dictionary of lists of phi/psi-angles for each residues
def get_phi_psi_angles(chain):

    # Use method from biopython to extract phi/psi angles
    polypeptide = Polypeptide.Polypeptide(chain)
    phi_psi_disordered_list = polypeptide.get_phi_psi_list()

    phi_psi_dictionary = dict()

    # Put everything into a nicely formatted dictionary
    for res_index, residue in enumerate(polypeptide):
        phi_psi_tuple = phi_psi_disordered_list[res_index]

        phi_psi_dictionary[residue.id[1]] = [phi_psi_tuple[0],
                                             phi_psi_tuple[1]]

    return phi_psi_dictionary
    

# Return a dictionary of lists of chi-angles for each residues
def get_chi_angles(chain):

    chi_dictionary = dict()

    # Put everything into a nicely formatted dictionary
    for residue in chain:
        chi_dictionary[residue.id[1]] = calculate_chi_angles(residue)

    return chi_dictionary
	

# Load first chain of first model found in the pdb-file
def load_chain(filename):

    for model in PDBParser().get_structure("chain", filename) :
        for chain in model:

            return chain




