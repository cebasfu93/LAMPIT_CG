import numpy as np
from optparse import OptionParser
import math, random
from  scipy.spatial.distance import cdist
from  transformations import *
from sklearn.decomposition import PCA
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

parser=OptionParser()
parser.add_option("-c", "--core", action="store", type='string', dest="CoreFile", default='NP1.xyz', help="Name of the core input file (xyz)")
parser.add_option("-e", "--element", action="store", type='string', dest="CoreElement", default='Pt', help="Chemical symbol of the metal in the core")
parser.add_option("-l", "--ligand", action="store", type='string', dest="LigandFile", default='LF1.gro', help="Name of the ligand input file (gro)")
parser.add_option("-a", "--anchor", action="store", type='int', dest="AnchorIndex", default='-1', help="Index (starting on 1) in the ligand's gro file corresponding to the atom connecting the core")
parser.add_option("-o", "--output", action="store", type='string', dest="OutputFile", default='NP1.pdb', help="Name of the output file (gro)")
(options, args) = parser.parse_args()
core_opt = options.CoreFile
element_opt = options.CoreElement
ligand_opt = options.LigandFile
anchor_opt = options.AnchorIndex-1
output_opt = options.OutputFile

def load_lig():
    lig_file = np.genfromtxt(ligand_opt, delimiter='\n', dtype='str', skip_header=2, skip_footer=1)
    xyz = []
    names = []
    resids = []
    for line in lig_file:
        xyz.append(str(line).split()[3:6])
        names.append(str(line).split()[1])
        resids.append(str(line).split()[0][-3:])
    xyz, names, resids = np.array(xyz, dtype='float'), np.array(names), np.array(resids)

    anchor_pos = xyz[anchor_opt]
    for i in range(len(xyz)):
        xyz[i,:] = xyz[i,:] - anchor_pos
    return xyz, names, resids

def load_core():
    core_file = np.genfromtxt(core_opt, delimiter='\n', dtype='str', skip_header=2)
    xyz = []
    names = []
    for line in core_file:
        xyz.append(str(line).split()[1:4])
        names.append(str(line).split()[0])
    xyz, names = np.array(xyz, dtype='float'), np.array(names)

    COM=np.average(xyz[names==element_opt,:], axis=0)
    for i in range(len(xyz)):
        xyz[i,:]=xyz[i,:]-COM
    return xyz/10., names

def get_ligand_pill(xyz_lig_func, anchor_ndx_func):
    #Runs a PCA and takes the first eigenvector as the best fitting line.
    pca = PCA(n_components=3)
    pca.fit(xyz_lig_func)
    pca1 = pca.components_[0]
    var1 = pca.explained_variance_[0]/np.sum(pca.explained_variance_)*100

    print("PCA1 explains: {:.1f}% of the points' variance".format(var1))

    #Randomly takes 2 other atoms in the ligand and project their positions in PCA1
    random.seed(666)
    rango = list(range(len(xyz_lig_func[:,0])))
    rango.remove(anchor_ndx_func)
    pillars_ndx = random.sample(rango, 2)
    pillars_func = np.array([0.0, 0.0, 0.0]) #This corresponds to the first stone (i.e. the anchor) which will always be in (0,0,0)
    for i in pillars_ndx:
        pillars_func = np.vstack((pillars_func, np.dot(xyz_lig_func[i], pca1) * pca1))
    return pillars_func

def get_stones(xyz_pillars):
    xyz_S_core = x_core[n_core == 'S']
    stones = []
    for i in range(len(xyz_S_core)):
        anchor_act = []
        norma = np.linalg.norm(xyz_S_core[i])
        for j in range(len(xyz_pillars)):
            D_S = np.linalg.norm(xyz_pillars[j])
            new_S = (norma + D_S)/norma * xyz_S_core[i]
            anchor_act.append(new_S)
        stones.append(anchor_act)
    return np.array(stones)

def coat_NP(names_lig, xyz_lig, xyz_stones, xyz_pillars):
    #Merges xyz coordinates and names of the core and the ligands into one coated NP
    xyz_coated = x_core[n_core==element_opt]
    names_coated = n_core[n_core==element_opt]

    #Transforms and appends rototranslated ligand 1
    xyz_lig_conv = np.insert(xyz_lig, 3, 1, axis=1).T
    for i in range(len(xyz_stones[:,0,0])):
        xyz_stones_now = xyz_stones[i,:,:]
        trans_matrix = affine_matrix_from_points(xyz_pillars.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig = np.dot(trans_matrix, xyz_lig_conv).T[:,:3]

        xyz_coated = np.append(xyz_coated, trans_lig, axis=0)
        names_coated = np.append(names_coated, names_lig, axis=0)
    return xyz_coated, names_coated

def print_gro(xyz_coat, names_coat, resid_coat):

    at = 0
    res = 0

    output=open(output_opt, "w")
    output.write(core_opt + " - " + ligand_opt + "\n")
    output.write(str(len(xyz_coat)) + "\n")
    output.close()
    for i in range(len(xyz_coat)):
        at += 1
        rname = resid_coat[i]
        if rname == element_opt:
            res += 1
        elif rname != resid_coat[i-1] or names_coat[i] == names_coat[i-1]:
            res += 1
        write_gro_block(res ,resid_coat[i], names_coat[i], at, xyz_coat[i], output_opt)
    output=open(output_opt, "a")
    output.write("\t 10.0 10.0 10.0\n")
    output.close()

def write_gro_block(resnum, resname, atname, at, xyz, ofile):
    output=open(ofile, "a")
    output.write('{:5d}{:5}{:>5}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(resnum, resname, atname, at, xyz[0], xyz[1], xyz[2]))
    output.close()

x_lig, n_lig, r_lig = load_lig()
x_core, n_core = load_core()
x_pillars = get_ligand_pill(x_lig, anchor_opt)
x_stones = get_stones(x_pillars)
x_coat, n_coat = coat_NP(n_lig, x_lig, x_stones, x_pillars)

r_coat = n_core[n_core==element_opt]
N_S = len(n_core[n_core=='S'])
for i in range(N_S):
    r_coat = np.append(r_coat, r_lig)

print_gro(x_coat, n_coat, r_coat)
#print_gro(x_core, n_core, n_core)
