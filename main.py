# Copyright (c) 2012 Anders S. Christensen
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

import structure_tools as struct
import sys

def calculate_chemical_shifts(chain):
    
    # Calculate all phi/psi angles once
    phi_psi_angles = struct.get_phi_psi_angles(chain)
    
    
    # Calculate all chi-angles once
    chi_angles = struct.get_chi_angles(chain)
    
    
    # Loop over residue in chain:
    for residue in chain:
    
        res_id = residue.id[1]
        res_name = residue.get_resname()
    
        # Do something ...
        #
        # print res_name, res_id, phi_psi_angles[res_id], chi_angles[res_id]


    # Retur something ...
    return []



def main():

    try:
        pdb_filename = sys.argv[1]
        
        # Load chain (as Bio.PDB chain-type) from pdb-file
        chain = struct.load_chain(pdb_filename)
        print 'Loaded chain from file:', pdb_filename

    except:
        print 'ERROR: Problem loading chain from file:', pdb_filename
        raise SystemExit


    # Calculate chemical shifts
    chemical_shifts = calculate_chemical_shifts(chain)




   
