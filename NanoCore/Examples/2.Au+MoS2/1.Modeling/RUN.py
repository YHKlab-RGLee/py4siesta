from NanoCore import *

# read structures from external files
mos2 = io.read_xsf('2H.xsf')    # 2d pbc on xz plane
slab = io.read_xsf('Au111.xsf') # 3d pbc 

# get cell matching
n1, n2 = slab.get_cell_match(mos2, 'y', 'y', max_unit=5)

# set repeat unit (USER INPUT)
m1, m2 = 3, 4

# generate supercells
slab_2 = slab * (m1, n1, 1)
mos2_2 = mos2 * (1, n2, m2)

# adjust cell size
ratio = mos2_2.get_cell()[1][1] / slab_2.get_cell()[1][1]
slab_3 = slab_2.adjust_cell_size(ratio, 7)

# set distance
d = 2.66328                    # arb. distance
zmax = slab_3.get_zmax()       # location of the surface
mos2_2.select_all()
mos2_2.translate(0, 0, zmax+d) # move MoS2 on the surface

# add two structures
new = slab_3 + mos2_2  # merge two systems with PBC of "slab"

# set vacuum along z axis
new.set_vacuum(35.0, 'z')

# sort atoms by z coordinates
new.select_all()
new.sort('z')

# save the interface model
io.write_xsf('model.xsf', new)
