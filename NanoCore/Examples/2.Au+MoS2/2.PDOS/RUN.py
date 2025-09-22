from NanoCore import *

# modeling tool
atoms = io.read_xsf('model.xsf')

# simulation object
sim = s2.Siesta(atoms)

# set simulation options
sim.set_option('kgrid', [1,9,1])
sim.set_option('kshift', [0.0,0.0,0.0])
sim.set_option('MixingWt', 0.02)
sim.set_option('PDOS', 1)
sim.set_option('PDOSE', (-5,5,0.1,1001))

# run siesta
sim.run(log=0)

# get pdos
E1, dos11, dos12 = s2.get_pdos(sim, -5, 5, by_atom=0, 
                               species=['Mo'], broad=0.05, npoints=1001)
E2, dos21, dos22 = s2.get_pdos(sim, -5, 5, by_atom=0, 
                               species=['S'], broad=0.05, npoints=1001)

# customized figure 
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(10,5))
fig1 = fig.add_subplot(111)

# draw dos
fig1.plot(E1, dos11, label='Mo')
fig1.plot(E2, dos21, label='S')
fig1.axis([-5, 5, 0, 12])
fig1.legend()
plt.show()
