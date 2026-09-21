"""
def Murnaghan(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    equation From PRB 28,5480 (1983)
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]

    E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)

    return E


def objective(pars,y,x):
    err =  y - Murnaghan(pars,x)
    return err
"""

def get_eos(pattern='*', struct_file='STRUCT.fdf'):
    # should be replaced by something using xml parser...
    dirs = glob(pattern)
    dirs.sort()

    volume = []
    for f in dirs:
        os.chdir(f)
        atoms = read_fdf("%s" % struct_file)
        cell = atoms.get_cell()
        v = Vector(cell[0]).dot( Vector(cell[1]).cross(Vector(cell[2])) )
        volume.append(v)
        os.chdir('..')
    os.system("grep 'siesta:         Total =' */stdout.txt > OUT")
    lines = open('OUT').readlines()
    volume = np.array(volume)

    energy = []
    for line in lines:
        e = float(line.split()[-1])
        energy.append(e)
    energy = np.array(energy)

    import pylab as plb # this includes numpy as np!
    from scipy.optimize import leastsq

    # make a vector to evaluate fits on with a lot of points so it looks smooth                    
    vfit = np.linspace(min(volume),max(volume),100)

    ### fit a parabola to the data
    # y = ax^2 + bx + c
    a,b,c = plb.polyfit(volume, energy, 2) #this is from pylab

    # now here are our initial guesses.
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4

    # now we have to create the equation of state function
    def Murnaghan(parameters,vol):
        '''
        given a vector of parameters and volumes, return a vector of energies.
        equation From PRB 28,5480 (1983)
        '''
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]
        E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
        return E

    # and we define an objective function that will be minimized
    def objective(pars,y,x):
        # we will minimize this function
        err =  y - Murnaghan(pars,x)
        return err

    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

    murnpars, ier = leastsq(objective, x0, args=(energy, volume)) #this is from scipy

    # now we make a figure summarizing the results
    plb.plot(volume,energy,'ro')
    plb.plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
    plb.plot(vfit, Murnaghan(murnpars,vfit), label='Murnaghan fit')
    plb.xlabel('Volume ($\AA^3$)')
    plb.ylabel('Energy (eV)')
    plb.legend(loc='best')

    # add some text to the figure in figure coordinates
    ax = plb.gca()
    plb.text(0.4,0.5,'Min volume = %1.2f $\AA^3$' % murnpars[3],
             transform = ax.transAxes)
    plb.text(0.4,0.4,'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (murnpars[1],
                                                                      murnpars[1]*160.21773)
             , transform = ax.transAxes)
    plb.savefig('a-eos.png')
    plb.show()

    print 'initial guesses  : ',x0
    print 'fitted parameters: ', murnpars

