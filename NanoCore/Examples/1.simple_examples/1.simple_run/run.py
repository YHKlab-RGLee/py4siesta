from NanoCore import *

n = 10
m = 0


# modeling tool
atom = carbonlab.cnt(n,m)
#atom = carbonlab.grp(0,0) # graphene
#atom = carbonlab.grp_nr(10,0) # graphene nanoribbons (zigzag)
#atom = carbonlab.cnt(n,m)


# Traslate the structrure to center of cell
atom.select_all()
center = atom.center(mode="geom")
cell = np.array(atom.get_cell())
distance = (cell[0]+cell[1]+cell[2])/2 - np.array(center)
atom.translate(*distance)

# simulation object
sim = s2.Siesta(atom)

# set simulation options
sim.set_option('kgrid', [1,1,100])      # set #kpoints
sim.set_option('kshift', [0,0,0.0]) # set k shift from gamma
sim.set_option('MixingWt', 0.10)        # adjust mixing weight (density)
sim.set_option('BasisSize', 'DZP')       # adjust basis size

# run siesta: 1st run
sim.run()
e1 = s2.get_total_energy()
print(e1)
