# Copyright (c) 2012, Anders S. Christensen
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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




   
