import numpy as np
from optparse import OptionParser
from  scipy.spatial.distance import cdist
import itertools


parser=OptionParser()
parser.add_option("-f", "--np", action="store", type='string', dest="NPFile", default='NP1.gro', help="Name of the coated nanoparticle input file (gro)")
parser.add_option("-o", "--output", action="store", type='string', dest="OutFile", default='NP1.top', help="Output file name")
(options, args) = parser.parse_args()
NP_opt = options.NPFile
out_opt = options.OutFile
M_bead = "C0"
cons_type = "1" #as written by martinize
bond_type = "1"
EN_cons = 35000 #32500

def init_gro(name_file):
    gro_file = np.genfromtxt(name_file, dtype='str', delimiter="\n", skip_header=2, skip_footer=1)
    xyz = []
    names = []
    resids = []
    for i in range(len(gro_file)):
        line = gro_file[i]
        xyz.append(line[-24:].split())
        names.append(line[-34:-29].split()[0])
        resids.append(line[-39:-34].split()[0])
    return np.array(xyz).astype('float'), np.array(names), np.array(resids)

def write_headers():
    out.write("; Topology written by Top_builder.py \n \n")
    out.write("#include \"inputs/martini_v2.2refP_SFU.itp\" \n \n")

def write_moltype():
    out.write("[ moleculetype ]\nCORE \t 3\n")

def write_atoms():
    out.write("[ atoms ] \n;   nr    type   resnr  residu    atom    cgnr  charge \n")
    for i in range(NM):
        out.write("{:5d}{:>6}{:6d}{:>6}{:>6}{:6d}{:8.4f} ; None \n".format(i+1, M_bead, i+1, "Pt", "Pt", i+1, 0.0))
    out.write("\n")

def write_bonds():
    out.write("[ bonds ] \n;  ai    aj funct           c0           c1 \n")
    #write_concentric_M_bonds()
    write_M_bonds()
    out.write("\n")

def write_concentric_M_bonds():
    D_M_M = cdist(x_M, x_M)
    norm_M = np.linalg.norm(x_M, axis=0)
    N_M = len(x_M)
    for i in range(1, N_M):
        out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; M - M \n".format(1, i+1, bond_type, D_M_M[0,i], EN_cons))

def write_M_bonds():
    D_M_M = cdist(x_M, x_M)
    bonds = []
    for i in range(NM):
        near_M = np.append(np.argsort(D_M_M[i])[1:7], np.argsort(D_M_M[i])[-1])
        #near_M = np.argsort(D_M_M[i])[1:3]
        for j in range(len(near_M)):
            at1 = i+1
            at2 = near_M[j]+1
            a = [at1, at2]
            if a not in bonds :
                a.reverse()
                bonds.append(a)
                out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; M - M \n".format(at1, at2, bond_type, D_M_M[at1-1, at2-1], EN_cons))

def write_footers():
    out.write("[ system ] \n; name \n" + str("CORE") + "\n \n")
    out.write("[ molecules ] \n; name \t \t number \n"+str("CORE") + "\t \t 1 \n")

x_M, n_M, r_M = init_gro(NP_opt)
NM = len(n_M)

out = open(out_opt,'w')

write_headers()
write_moltype()
write_atoms()
write_bonds()
write_footers()

out.close()
print("DONT FORGET TO CHANGE THE MASS OF C0 IN THE MARTINI_SFU ITP FILE!")
