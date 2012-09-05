import sys
import structure_tools as struct


# Get input pdb-file
pdb_filename = sys.argv[1]


# Load chain (as Bio.PDB chain-type) from pdb-file
chain = struct.load_chain(pdb_filename)


# Calculate all phi/psi angles once
phi_psi_angles = struct.get_phi_psi_angles(chain)


# Calculate all chi-angles once
chi_angles = struct.get_chi_angles(chain)


# Loop over residue in chain:
for residue in chain:

    res_id = residue.id[1]
    res_name = residue.get_resname()

    print res_name, res_id, phi_psi_angles[res_id], chi_angles[res_id]

    # for atom in residue:
        # print atom

