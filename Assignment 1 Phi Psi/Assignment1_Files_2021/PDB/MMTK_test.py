#!/usr/bin/python

"""
Structural Bioinformatics assignment 1 - Phi-psi angles

When you finish the assignment, the script should create an
output file containing the phi and psi angles and secondary
structure assignment of each residue. Please ONLY modify the
code in the three indicated blocks and do NOT use additional
python packages. Use spaces instead of tabs. Do not round
down in any calculation step. The commented #print() lines
can be used to test separate functions, but make sure that
all of them are commented when submtting to CodeGrade.

To run, make 'PDB' your working directory and use:

> python3 readPDB.py pdb_filename.txt
"""

# Packages
from sys import argv
import os
from math import sqrt, atan2, degrees

# Vector functions that we need to calculate the angles
def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
#print(dot_product([1, 2, 3], [1, 3, 2]))

def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1]*v2[2] - v1[2]*v2[1]
    j = v1[2]*v2[0] - v1[0]*v2[2]
    k = v1[0]*v2[1] - v1[1]*v2[0]
    return [i,j,k]
#print(cross_product([1, 2, 3], [1, 3, 2]))

def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)
#print(magnitude([1, 2, 2]))

# PDB file parser
def readPDB(PDB_file):
    """ Reads a PDB file and stores the atom
    coordinates and amino acid types of the protein """
    # open the file
    f = open(PDB_file, 'r')

    # dictionaries to store the output
    # pdb atom coordinates:
    #     pdbcoord[chain][residue_number][atom_type] = coordinates
    pdbcoord = {}
    # residue type per chain and residue number (i.e. store the sequence)
    #     pdbseq[chain][resnum] = restype
    pdbseq = {}

    # parse each line in the file
    for line in f:
        # remove whitespace at the end of the line
        line = line.strip()
        # only parse the lines containing atom coordinates
        if line[:4] == 'ATOM':
            # ATOM type (e.g. C-alpha)
            atom_type = line[12:16].strip()
            # AMINO ACID type (e.g. alanine)
            aa_type = line[17:20].strip()
            # residue number
            res_num = int(line[22:26])
            # Protein chain
            chain = line[21]
            # coordinates
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])

            # if chain does not exists create new entry
            if not chain in pdbcoord:
                pdbcoord[chain] = {}
                pdbseq[chain] = {}
            # if resnum does not exists create new entry
            if not res_num in pdbcoord[chain]:
                pdbcoord[chain][res_num] = {}

            # store coordinates as a vector
            pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]
            # store sequence
            pdbseq[chain][res_num] = aa_type

    # close file
    f.close()

    # return dictionaries
    return pdbcoord, pdbseq
#print(readPDB('Data\\1TIM.pdb'))

### THE FOLLOWING THREE FUNCTIONS ARE THE ONES YOU NEED
### TO EDIT FOR THE ASSIGNMENT. ONLY EDIT THE INDICATED
### BLOCKS
def calculateDihedral(a1, a2, a3, a4):
    """ Calculates the normal vector of the planes
    defined by four atom coordinates """

    ### START CODING HERE
    # calculate normal vectors to the planes defined by a1,a2,a3 and a2,a3,a4
    # you may use the functions "cross_product","dot_product" and "magnitude" defined above
    # you can also use the python math function "atan2" and "degrees"

    # calculate connecting vectors (bonds)
    #Get vectors b1 (=a2-a1), b2 (=a3-a2) and b3 (=a4-a3)
    b1 = [a2[0]-a1[0],a2[1]-a1[1],a2[2]-a1[2]]
    b2 = [a3[0]-a2[0],a3[1]-a2[1],a3[2]-a2[2]]
    b3 = [a4[0]-a3[0],a4[1]-a3[1],a4[2]-a3[2]]
    print(b1, b2, b3)

    n1 = [3,4,0] #cross_product(b1,b2)
    n2 = [1,2,3] #cross_product(b2,b3)
    print(n1, n2)

	#u1
    u1 = [b2[0]/magnitude(b2), b2[1]/magnitude(b2), b2[2]/magnitude(b2) ]
    print(u1)

    #cos_dihedral
    dot = dot_product(n1,n2)
    mag = (magnitude(n1)*magnitude(n2))
    cos_dihedral =  dot / mag
    print(dot, mag, cos_dihedral)

    #sin_dihedral
    sin_top_1 = magnitude(b2)*b1[0] + magnitude(b2)*b1[1] + magnitude(b2)*b1[2]
    sin_top_2 = dot_product(sin_top_1, n2)
    sin_dihedral = sin_top_2 / mag
    sin_dihedral = cos_dihedral*u1[0] + cos_dihedral*u1[1] + cos_dihedral*u1[2]
    print("sin", sin_dihedral)
    dihedral = degrees(atan2(sin_dihedral, cos_dihedral))

    #dihedral = 0 # replace this line with your code
    ### END CODING HERE
    return dihedral
print(calculateDihedral([1, 9, 2], [3, 2, 1], [2, 4, 7], [8, 2, 5]))

def assign_ss(phi, psi):
    """ Assign a secondary structure type based on the phi
    and psi angles of a residue """
    ### START CODING HERE
    # for code checking purposes use the terms "loop", "alpha" or "beta"
    #secondary_structure = "" # replace this line with your code
    if (90 <= phi and phi<= 180 and -180 <= psi and psi <= -50):
        secondary_structure = "beta"
    elif (-60 <= phi and phi <= -40 and -160 <= psi and psi<= -50):
        secondary_structure = "alpha"
    else:
        secondary_structure = "loop"
    ### END CODING HERE
    return secondary_structure
#print(assign_ss(60, 25))

def print_phi_psi(pdbcoord, pdbseq, outfile):
    """ given the PDB coordinates, calculate the dihedral
    angles of all the residues, assign secondary structure
    types and write them into an output file """
    f = open(outfile, 'w')

    # get the chains from the PDB file
    list_chains = sorted(pdbcoord.keys())

    for chain in list_chains:
        # get the sorted residue numbers from the pdbcoord dictionary
        list_residue_numbers = sorted(pdbcoord[chain].keys())
        for res_num in list_residue_numbers:
            # if certain residues are missing in the PDB file, you will
            # get a KeyError. Make sure your program does not crash, but
            # gives a warning when this happens
            try:
                ### START CODING HERE
                phi = psi = ss = None
    			### START CODING HERE (and calculate phi, psi and ss)
    			#Get coordinates for N, Ca an CO of the res_num
                N = pdbcoord[chain][res_num]["N"]
                CA = pdbcoord[chain][res_num]["CA"]
                CO = pdbcoord[chain][res_num]["C"]
    			#Get coordinates for C_prev of previous res_num
                CO_prev = pdbcoord[chain][res_num-1]["C"]
    			#Now C_prev, N, CA and C can be used to calc phi for res_num
                #print(CO_prev,N,CA,CO)
                phi = 0 #calculateDihedral(CO_prev,N,CA,CO)
                #print(phi)

    			#Find atom_num for N_next of next res_num
                N_next = pdbcoord[chain][res_num+1]["N"]
                #Now N, CA, CO and N_next can be used to calc psi for res_num
                psi = 0 #calculateDihedral(N,CA,CO,N_next)
    			#print(phi,psi)

                #Assign secondary structure based on phi and PSI
                ss = assign_ss(phi,psi)

                ### END CODING HERE
                #
                #phi, psi, ss = 0, 0, "test" # replace this line with your code

            except KeyError:
                print('WARNING: KeyError:', KeyError, 'in residue', chain, res_num)

            # get amino acid
            aa_type = pdbseq[chain][res_num]
            # write into output file
            print(chain, res_num, aa_type, phi, psi, ss, file=f)
    f.close()
    print('written:', outfile)

def main():
    # input PDB file
    f_in = argv[1]
    f_out = 'Output/phi_psi.txt'
    # read PDB file
    pdbcoord, pdbseq = readPDB(f_in)
    print_phi_psi(pdbcoord, pdbseq, f_out)
    # for testing
    # for i in ['1TIM', '3PG8']:
    #     f_in = 'student/{}.pdb'.format(i)
    #     print(f_in)
    #     f_out = 'student/output/phi_psi_{}.txt'.format(i)
    #
    #     # read PDB file
    #     pdbcoord, pdbseq = readPDB(f_in)
    #     print_phi_psi(pdbcoord, pdbseq, f_out)

if __name__ == '__main__':
    main()
