import yaml
import numpy as np

data = yaml.load(open("qpoints.yaml"))
dynmat = []
dynmat_data = data['phonon'][0]['dynamical_matrix']
for row in dynmat_data:
    vals = np.reshape(row, (-1, 2))
    dynmat.append(vals[:, 0] + vals[:, 1] * 1j)
dynmat = np.array(dynmat)
print (dynmat.shape)

eigvals, eigvecs, = np.linalg.eigh(dynmat)
frequencies = np.sqrt(np.abs(eigvals.real)) * np.sign(eigvals.real)
conversion_factor_to_THz = 15.633302

for j in frequencies:
    print (j)
#    print (j*conversion_factor_to_THz)
#print (frequencies * conversion_factor_to_THz)
