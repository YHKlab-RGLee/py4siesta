from NanoCore import *

# modeling tool
atoms = io.read_xsf('model.xsf')

# simulation object
sim = s2.Siesta(atoms)

# set simulation options
sim.set_option('kgrid', [1,7,1])
sim.set_option('kshift', [0.0,0.0,0.0])
sim.set_option('MixingWt', 0.05)
sim.set_option('VH', 1)

# run siesta
sim.run(log=0)

# post-process by siesta scripts
z, v = s2.get_hartree_pot_z()

# customized plot
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(10,5))
fig1 = fig.add_subplot(111)
fig1.plot(z, v)
plt.show()
