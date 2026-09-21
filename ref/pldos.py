from nanocore import *
from matplotlib import cm, colors
import time, sys,os,glob
import os

#How to use
#python pldos.py Device

fname = sys.argv[1] # python pldos.py system_name
at = io.read_xyz('%s.xyz' % fname) # % label)
atoms = at._atoms

# slice atoms by z coordinates
z_coords = []; indice = []
for atom in atoms:
    if not atom[2] in z_coords: z_coords.append(atom[2])

for z in z_coords:
    temp = []
    for atom in atoms:
        if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
    indice.append(temp)


#find fermilevel
a=glob.glob('%s.EIG' % fname)
a = str(a[0])
f=open(a)
list_lines=[]
for line in f.readlines():
    list_lines.append(line)
Fermi = list_lines[0]
Fermi = Fermi.split()
Fermi = float(Fermi[0])

print( Fermi )
# get pdos
Z = []; E = []
for ind in indice:
    E1, dos11, dos12 = siesta.get_pdos(0.0, -10.0, 0.0, by_atom=1, atom_index=ind, broad= 0.02, npoints = 1001, label = fname)
    E = np.array(E1)
    Z.append(np.array(dos11))

absZ = np.abs(Z).T


# convert to log10 values + minimum correction to avoid -INF
Z = np.log10(absZ + 10**-5)

#Z = absZ

# generate meshgrid
X, Y = np.meshgrid(np.array(z_coords), E-Fermi)

# customized figure
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.xaxis.set_tick_params(width=1.5)
ax.yaxis.set_tick_params(width=1.5)

ax.tick_params(axis="x", direction='in', length=12)
ax.tick_params(axis="y", direction='in', length=12)


Zmin = -4
Zmax = 2

plt.ylim(-4, 2)
plt.xlim(25, 85)

plt.yticks(np.arange(-4,2,1))
plt.xticks(np.linspace(25,85,5))

# labels
labels = [item.get_text() for item in ax.get_xticklabels()]
empty_string_labels = ['']*len(labels)
ax.set_xticklabels(empty_string_labels)
labels = [item.get_text() for item in ax.get_yticklabels()]
empty_string_labels = ['']*len(labels)
ax.set_yticklabels(empty_string_labels)



levels = np.linspace(Zmin, Zmax, 100)
print(Z.min(), Z.max())


cmap = plt.cm.get_cmap("viridis")


import pylab as plb
cset = plb.contourf(X, Y, Z, levels, cmap=cmap, extend = 'both')

# add colorbar
cbar = fig.colorbar(cset, ticks=np.linspace(Zmin,Zmax,3), extend = 'both')
cbar.ax.tick_params(labelsize=14)

plt.show()
fig.savefig('pldos.png', dpi=300, transparent = True)


