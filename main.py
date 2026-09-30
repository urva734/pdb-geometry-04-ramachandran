import math
import matplotlib.pyplot as plt

pdb_file = "1aki.pdb"
residues = {}
with open(pdb_file) as f:
    for line in f:
        if not line.startswith("ATOM"): continue
        atom = line[12:16].strip()
        if atom not in ("N","CA","C"): continue
        res_id = int(line[22:26].strip())
        x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
        if res_id not in residues: residues[res_id]={}
        residues[res_id][atom]=(x,y,z)

def sub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def mag(v): return math.sqrt(dot(v,v))
def norm(v):
    m=mag(v)
    return (v[0]/m, v[1]/m, v[2]/m) if m>0 else v

def calc_dihedral(p1,p2,p3,p4):
    b1=sub(p2,p1); b2=sub(p3,p2); b3=sub(p4,p3)
    n1=cross(b1,b2); n2=cross(b2,b3)
    n1=norm(n1); n2=norm(n2); b2n=norm(b2)
    m1=cross(n1,b2n)
    x=dot(n1,n2); y=dot(m1,n2)
    return math.degrees(math.atan2(y,x))

def classify(phi,psi):
    if -140 <= phi <= -20 and -80 <= psi <= 10:
        return "ALPHA"
    if -180 <= phi <= -40 and 90 <= psi <= 180:
        return "BETA"
    if -180 <= phi <= -40 and -180 <= psi <= -90:
        return "BETA"
    if -150 <= phi <= -20 and 20 <= psi <= 90:
        return "BETA"
    if 20 <= phi <= 120 and -30 <= psi <= 80:
        return "LEFT/GLY"
    return "OUTLIER"

print("=== PDB Geometry 04: Ramachandran Check ===")
print(f"File: {pdb_file} | Total residues: {len(residues)}\n")
ids=sorted(residues.keys())
alpha=beta=left=outlier=0
phi_a=[]; psi_a=[]; phi_b=[]; psi_b=[]; phi_l=[]; psi_l=[]; phi_o=[]; psi_o=[]

for i in range(1,len(ids)-1):
    prev_id=ids[i-1]; curr_id=ids[i]; next_id=ids[i+1]
    if 'C' not in residues.get(prev_id,{}): continue
    if not all(k in residues.get(curr_id,{} ) for k in ('N','CA','C')): continue
    if 'N' not in residues.get(next_id,{}): continue
    phi=calc_dihedral(residues[prev_id]['C'],residues[curr_id]['N'],residues[curr_id]['CA'],residues[curr_id]['C'])
    psi=calc_dihedral(residues[curr_id]['N'],residues[curr_id]['CA'],residues[curr_id]['C'],residues[next_id]['N'])
    region=classify(phi,psi)
    if region=="ALPHA":
        alpha+=1; tag="ALPHA [OK]"; phi_a.append(phi); psi_a.append(psi)
    elif region=="BETA":
        beta+=1; tag="BETA [OK]"; phi_b.append(phi); psi_b.append(psi)
    elif region=="LEFT/GLY":
        left+=1; tag="LEFT [OK]"; phi_l.append(phi); psi_l.append(psi)
    else:
        outlier+=1; tag="OUTLIER [VIOLATION]"; phi_o.append(phi); psi_o.append(psi)
    print(f"Residue {curr_id:3d} : Phi {phi:7.1f}, Psi {psi:7.1f} : {tag}")

print(f"\nTotal checked: {alpha+beta+left+outlier}")
print(f"Alpha: {alpha} | Beta: {beta} | Left/Gly: {left} | Outliers: {outlier}")

# IMAGE
plt.figure(figsize=(7,6))
plt.scatter(phi_a, psi_a, c='red', label=f'Alpha ({alpha})')
plt.scatter(phi_b, psi_b, c='blue', label=f'Beta ({beta})')
plt.scatter(phi_l, psi_l, c='green', label=f'Left ({left})')
plt.scatter(phi_o, psi_o, c='black', label=f'Outlier ({outlier})')
plt.xlim(-180,180); plt.ylim(-180,180)
plt.xlabel("Phi"); plt.ylabel("Psi")
plt.title("Ramachandran Plot - 1aki.pdb")
plt.legend(); plt.grid(alpha=0.3)
plt.savefig("ramachandran_plot.png", dpi=150)
print("Saved: ramachandran_plot.png")