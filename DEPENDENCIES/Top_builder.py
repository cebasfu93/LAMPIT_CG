import numpy as np
from optparse import OptionParser
from  scipy.spatial.distance import cdist
import itertools

parser=OptionParser()
parser.add_option("-f", "--np", action="store", type='string', dest="NPFile", default='NP1.gro', help="Name of the coated nanoparticle input file (gro)")
parser.add_option("-t", "--title", action="store", type='string', dest="NPName", default='NP1', help="Name of the coated nanoparticle input file")
parser.add_option("-l", "--lig", action="store", type='string', dest="LigFile", default='Lig.gro', help="Name of the ligand structure file (gro)")
parser.add_option("-i", "--itp", action="store", type='string', dest="ItpFile", default='Lig_A.itp', help="Name of the ligand's itp topology written by martinize")
parser.add_option("-p", "--top", action="store", type='string', dest="TopFile", default='Lig.top', help="Name of the ligand's top topology written by martinize")
parser.add_option("-e", "--element", action="store", type='string', dest="Element", default='Pt', help="Element")
parser.add_option("-a", "--anchor", action="store", type='int', dest="AnchorIndex", default='0', help="Index (starting on 1) in the ligand's gro file corresponding to the atom connecting the core")
parser.add_option("-o", "--output", action="store", type='string', dest="OutFile", default='NP1.top', help="Output file name")
(options, args) = parser.parse_args()
NP_opt = options.NPFile
NPTitle = options.NPName
Lig_opt = options.LigFile
Itp_opt = options.ItpFile
Top_opt = options.TopFile
Ele_opt = options.Element
anchor_opt = options.AnchorIndex-1
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

def init_topology(name_file):
    top_file = np.genfromtxt(name_file, dtype='str', delimiter="\n")
    return top_file

def write_headers():
    out.write("; Topology written by Top_builder.py \n \n")
    out.write("#include \"inputs/martini_v2.2refP_SFU.itp\" \n \n")
    out.write("#define RUBBER_BANDS \n \n")

def write_moltype():
    found_pattern = False
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ moleculetype ]" in line:
            found_pattern = True
        elif "[ atoms ]" in itp[i+1]:
            out.write(NPTitle + "\t \t 1 \n \n")
            break
        if found_pattern:
            out.write(line + "\n")

def write_atoms():
    found_pattern = False
    section = []
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ atoms ]" in itp[i]:
            found_pattern = True
        elif "[ virtual_sites2 ]" in itp[i]:
            break
        elif found_pattern and line[0] != ";":
            section.append(line.split())
    out.write("[ atoms ] \n;   nr    type   resnr  residu    atom    cgnr  charge \n")
    for i in range(NM):
        out.write("{:5d}{:>6}{:6d}{:>6}{:>6}{:6d}{:8.4f} ; None \n".format(i+1, M_bead, i+1, Ele_opt, Ele_opt, i+1, 0.0))
    #for i in range(6):
        #out.write("{:5d}{:>6}{:6d}{:>6}{:>6}{:6d}{:8.4f} ; None \n".format(NM+i+1, 'P5', NM+i+1, 'S', 'SC1', i+1, 0.0))

    section = np.array(section)
    res = NM
    resid = Ele_opt

    for i in range(NL):
        for j in range(len(section)):
            if resid != section[j][3]:
                res += 1
                resid = section[j][3]
            if len(section[j])==9:
                    out.write("{:5d}{:>6}{:6d}{:>6}{:>6}{:6d}{:8.4f} ; {} \n".format(NM+i*N_at_lig+j+1, section[j][1], res, resid, section[j][4], NM+i*N_at_lig+j+1, float(section[j][6]), section[j][8]))
            if len(section[j])==10:
                    out.write("{:5d}{:>6}{:6d}{:>6}{:>6}{:6d}{:8.4f}{:8.4f} ; {} \n".format(NM+i*N_at_lig+j+1, section[j][1], res, resid, section[j][4], NM+i*N_at_lig+j+1, float(section[j][6]), float(section[j][7]), section[j][9]))
    out.write("\n")

def write_virtual_sites2():
    found_pattern = False
    section = []
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ virtual_sites2 ]" in itp[i]:
            found_pattern = True
        elif "[ bonds ]" in itp[i]:
            break
        elif found_pattern and line[0] != ";":
            section.append(line.split())

    section = np.array(section)
    out.write("[ virtual_sites2 ] \n")
    for i in range(NL):
        for j in range(len(section)):
            at1 = NM+i*N_at_lig+int(section[j,0])
            at2 = NM+i*N_at_lig+int(section[j,1])
            at3 = NM+i*N_at_lig+int(section[j,2])
            out.write("{:5d}{:6d}{:6d}{:>7}{:10.5f} ; {} \n".format(at1, at2, at3, section[j,3], float(section[j,4]), section[j,6]))
    out.write("\n")

def write_bonds():
    found_pattern = False
    section = []
    section_rub = []
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ bonds ]" in itp[i]:
            found_pattern = True
        elif "[ constraints ]" in itp[i]:
            break
        elif found_pattern and line[0] != ";":
            if "RUBBER" not in line:
                section.append(line.split())
            else:
                section_rub.append(line.split())

    section, section_rub = np.array(section), np.array(section_rub)
    out.write("[ bonds ] \n;  ai    aj funct           c0           c1 \n")
    for i in range(NL):
        for j in range(len(section)):
            at1 = NM+i*N_at_lig+int(section[j,0])
            at2 = NM+i*N_at_lig+int(section[j,1])
            out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; {} \n".format(at1, at2, section[j,2], float(section[j,3]), int(section[j,4]), section[j,6]))

    #write_M_bonds()
    write_concentric_M_bonds()
    write_S_bonds()
    out.write("#ifndef NO_RUBBER_BANDS \n#ifndef RUBBER_FC \n#define RUBBER_FC 500.000000 \n#endif \n")
    for i in range(NL):
        for j in range(len(section_rub)):
            at1 = NM+i*N_at_lig+int(section_rub[j,0])
            at2 = NM+i*N_at_lig+int(section_rub[j,1])
            out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; RUBBER \n".format(at1, at2, section_rub[j,2], float(section_rub[j,3]), 500))
    out.write("#endif \n")
    out.write("\n")

def write_constraints():
    found_pattern = False
    section = []
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ constraints ]" in itp[i]:
            found_pattern = True
        elif "[ angles ]" in itp[i]:
            break
        elif found_pattern and line[0] != ";":
            section.append(line.split())

    section = np.array(section)
    out.write("[ constraints ] \n;  ai    aj funct           value \n")
    for i in range(NL):
        for j in range(len(section)):
            at1 = NM+i*N_at_lig+int(section[j,0])
            at2 = NM+i*N_at_lig+int(section[j,1])
            out.write("{:5d}{:6d}{:>7}{:10.5f} ; {} \n".format(at1, at2, section[j,2], float(section[j,3]), section[j,5]))

    #write_M_cons()
    #write_concentric_M_cons()
    #write_S_cons()
    out.write("\n")

def write_concentric_M_bonds():
    x_M = x_sys[n_sys==Ele_opt]
    D_M_M = cdist(x_M, x_M)
    norm_M = np.linalg.norm(x_M, axis=0)
    N_M = len(x_M)
    for i in range(1, N_M):
        out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; M - M \n".format(1, i+1, bond_type, D_M_M[0,i], EN_cons))

def write_M_bonds():
    D_M_M = cdist(x_sys[n_sys==Ele_opt], x_sys[n_sys==Ele_opt])
    bonds = []
    for i in range(NM):
        near_M = np.append(np.argsort(D_M_M[i])[1:5], np.argsort(D_M_M[i])[-1])
        #near_M = np.argsort(D_M_M[i])[1:5]
        for j in range(len(near_M)):
            at1 = i+1
            at2 = near_M[j]+1
            a = [at1, at2]
            if a not in bonds :
                a.reverse()
                bonds.append(a)
                out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; M - M \n".format(at1, at2, bond_type, D_M_M[at1-1, at2-1], EN_cons))

def write_M_cons():
    D_M_M = cdist(x_sys[n_sys==Ele_opt], x_sys[n_sys==Ele_opt])
    bonds = []
    for i in range(NM):
        near_M = np.argsort(D_M_M[i])[1:5]
        for j in range(len(near_M)):
            at1 = i+1
            at2 = near_M[j]+1
            a = [at1, at2]
            if a not in bonds :
                a.reverse()
                bonds.append(a)
                out.write("{:5d}{:6d}{:>7}{:10.6f} ; M - M \n".format(at1, at2, cons_type, D_M_M[at1-1, at2-1]))

def write_S_bonds():
    #D_S_M = cdist(x_sys[n_sys==n_lig[anchor_opt]], x_sys[n_sys==Ele_opt])
    for i in range(NL):
        D_S_M = cdist([x_sys[NM+i*N_at_lig+anchor_opt]], x_sys[n_sys==Ele_opt])[0]
        near_M = np.argsort(D_S_M)[0]
        at1 = NM+i*N_at_lig+anchor_opt+1
        at2 = near_M+1
        out.write("{:5d}{:6d}{:>7}{:10.5f}{:7d} ; M - S \n".format(at1, at2, bond_type, 0.47, EN_cons))

def write_S_cons():
    #D_S_M = cdist(x_sys[n_sys==n_lig[anchor_opt]], x_sys[n_sys==Ele_opt])
    for i in range(NL):
        D_S_M = cdist([x_sys[NM+i*N_at_lig+anchor_opt]], x_sys[n_sys==Ele_opt])[0]
        near_M = np.argsort(D_S_M)[0]
        at1 = NM+i*N_at_lig+anchor_opt+1
        at2 = near_M+1
        out.write("{:5d}{:6d}{:>7}{:10.5f} ; M - S \n".format(at1, at2, cons_type, 0.47))

def write_angles():
    found_pattern = False
    section = []
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ angles ]" in itp[i]:
            found_pattern = True
        elif "[ dihedrals ]" in itp[i]:
            break
        elif found_pattern and line[0] != ";":
            section.append(line.split())

    #section = np.array(section)
    out.write("[ angles ]\n;  ai    aj    ak funct           c0           c1\n")
    for i in range(NL):
        for j in range(len(section)):
            at1 = NM+i*N_at_lig+int(section[j][0])
            at2 = NM+i*N_at_lig+int(section[j][1])
            at3 = NM+i*N_at_lig+int(section[j][2])
            out.write("{:5d}{:6d}{:6d}{:>7}{:11.5f}{:6d} ; {} \n".format(at1, at2, at3, section[j][3], float(section[j][4]), int(section[j][5]), section[j][7]))
    out.write("\n")

def write_dihedrals():
    found_pattern = False
    section = []
    for i in range(len(itp)):
        line = str(itp[i])
        if "[ dihedrals ]" in itp[i]:
            found_pattern = True
        elif found_pattern == True and line[0] != ";":
            section.append(line.split())

    section = np.array(section)
    out.write("[ dihedrals ]\n ;  ai    aj    ak    al funct           c0           c1           c2 \n")
    for i in range(NL):
        for j in range(len(section)):
            at1 = NM+i*N_at_lig+int(section[j][0])
            at2 = NM+i*N_at_lig+int(section[j][1])
            at3 = NM+i*N_at_lig+int(section[j][2])
            at4 = NM+i*N_at_lig+int(section[j][3])
            if len(section[j])==10:
                out.write("{:5d}{:6d}{:6d}{:6d}{:>7}{:7d}{:6d}{:>7} ; {} \n".format(at1, at2, at3, at4, section[j][4], int(section[j][5]), int(section[j][6]), section[j][7], section[j][9]))
            if len(section[j])==9:
                out.write("{:5d}{:6d}{:6d}{:6d}{:>7}{:7d}{:6d} ; {} \n".format(at1, at2, at3, at4, section[j][4], int(section[j][5]), int(section[j][6]), section[j][8]))
    out.write("\n")

def write_footers():
    out.write("[ system ] \n; name \n" + str(NPTitle) + "\n \n")
    out.write("[ molecules ] \n; name \t \t number \n"+str(NPTitle) + "\t \t 1 \n")

x_sys, n_sys, r_sys = init_gro(NP_opt)
x_lig, n_lig, r_lig = init_gro(Lig_opt)
top = init_topology(Top_opt)
itp = init_topology(Itp_opt)
NM = len(n_sys[n_sys==Ele_opt])
N_at_lig = len(x_lig)
NL = int((len(x_sys)-NM)/N_at_lig)

out = open(out_opt,'w')

write_headers()
write_moltype()
write_atoms()
write_virtual_sites2()
write_bonds()
write_constraints()
write_angles()
write_dihedrals()
write_footers()

out.close()
print("DONT FORGET TO CHANGE THE MASS OF C0 IN THE MARTINI_SFU ITP FILE!")
