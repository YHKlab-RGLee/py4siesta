from NanoCore import *

# modeling tool
atoms = carbonlab.grp(0,0)

# simulation object
sim = s2.Siesta(atoms)

# set simulation options
sim.set_option('kgrid', [100,100,1])    # set #kpoints
sim.set_option('kshift', [0.5,0.5,0.0]) # set k shift from gamme
sim.set_option('MixingWt', 0.10)        # adjust mixing weight (density)
sim.set_option('BasisSize', 'SZ')       # adjust basis size

# run siesta
sim.run()

# get dos
emin = -20; emax = 5 # energy range for band & dos
sim.set_option('kgrid', [200,200,1])    # increase #kpoints
sim.run()                               # RE-run
E, dos, dos1, dos2 = s2.get_dos(emin=emin, emax=emax, broad=0.1)

# get reciprocal lattice vectors
v1, v2, v3 = atoms.get_reciprocal_cell(unit=np.pi)

# calculate special kpoints
g = np.array([0,0,0])
m = (v1 + v2)*0.5
k = (2*v1 + v2)*1./3.

# measure the length between each two special kpoints
l_kg = Vector(k-g).length()
l_gm = Vector(g-m).length()
l_mk = Vector(m-k).length()

# write bandpath into a file
nkpt = 100 # 1 kpoint / ang**-1
f_band = open('bandline', 'w')
f_band.write('BandLinesScale    pi/a\n')
f_band.write('%block BandLines\n')

for kp, lp in [(k, 1), (g, int(l_kg*nkpt)), (m, int(l_gm*nkpt)), (k, int(l_mk*nkpt))]:
    line = '%4i %8.4f %8.4f %8.4f\n' % (lp, kp[0], kp[1], kp[2])
    print (line)
    f_band.write(line)
f_band.write('%endblock BandLines\n')
f_band.close()

# get bandlines
path, eigs = s2.get_band(sim, './bandline')

# customized figure 
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(10,5))
fig1 = fig.add_subplot(121)
fig2 = fig.add_subplot(122)

# draw dos
fig2.plot(dos, E)
fig2.axis([0, np.array(dos).max(), emin, emax])

# draw band 
i = 0
for eig in eigs:
    fig1.plot(path[i], eig)
    i += 1
fig1.axis([0, np.array(path).max(), emin, emax])

plt.show()
