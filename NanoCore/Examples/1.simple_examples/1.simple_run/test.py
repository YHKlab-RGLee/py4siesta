from NanoCore import *

# Write Structure
atom = carbonlab.grp(0,0) # graphene
#atom = carbonlab.grp_nr(10,0) # graphene nanoribbons (zigzag)
#atom = carbonlab.cnt(10,0)

# Traslate the structrure to center of cell
atom.select_all()
center = atom.center(mode="geom")
cell = np.array(atom.get_cell())
distance = (cell[0]+cell[1]+cell[2])/2 - np.array(center)
atom.translate(*distance)



sys = s2.Siesta(atom)
sys.write_struct()
