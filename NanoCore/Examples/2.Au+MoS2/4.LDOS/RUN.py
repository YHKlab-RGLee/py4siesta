from NanoCore import *

# modeling tool
atoms = io.read_xsf('model.xsf')

# simulation object
sim = s2.Siesta(atoms)

# set simulation options
sim.set_option('kgrid', [1,9,1])
sim.set_option('kshift', [0.0,0.0,0.0])
sim.set_option('MixingWt', 0.05)
sim.set_option('LDOS', 1)            # activate LDOS option
sim.set_option('LDOSE', (-1.0, 1.0)) # LDOS energy range (min, max)

# run siesta
sim.run()

# post-process by siesta scripts
v1, v2, v3 = atoms.get_cell()    # define grid range
v2 = 3*Vector(v2)

origin = Vector(-0.2,-0.2,-0.2)  # define the origin
origin -= 0.5*Vector(v1)

nmesh = [150,150,350]             # the number of gridpoints along v1/2/3

s2.get_ldos(v1, v2, v3, origin, nmesh)
