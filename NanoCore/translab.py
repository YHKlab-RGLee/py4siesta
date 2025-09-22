from __future__ import print_function
from atoms import *
import siesta2 as s2


class Transport(object):
    """
Transport(elec1, chan, elec2=None)
    
    Class for representing a quantum transport simulation

    Parameters
    ----------
    electrode1/2 : AtomsSystem
        an atomic structure whose PBC is imposed along the transport direction
    channel : AtomsSystem 
        an atomic structure 

    Optional parameters
    -------------------

    Examples
    --------
    >>> elec = carbonlab.grp_armchair(4,4)
    >>> chan = carbonlab.chain(10, 1.285, 0.1)
    >>> trans_obj = Transport(elec, chan)
    
    """

    __slots__ = ['_elec1', '_elec2', '_chan']

    #
    # Initialization
    #
    def __init__(self, elec1, chan, elec2=None,
                 det_k=0, max_nk=10, init_nk=1, step_nk=1, opts=None):

        # 0. set the location of structures at z=0??

        # 1. set electrode and channel
        self.set_electrode1(elec1)
        if not elec2: self._elec2 = self._elec1.copy()
        self.set_electrode2(elec2)
        self.set_channel(chan)

        # 2. check: cell shape & size
        possible = self.is_possible()
        #if not possible: raise ValueError, "Please check your initial structures."

        # 3. determine number of kpoint
        if det_k: 
            kpts, energies = self.determine_nkpt(self._elec1,
                                                 max_nk, init_nk, step_nk, opts)
            i = 0
            for kpt in kpts:
                print (kpt, energies[i])
                i += 1

        # 4. make junction model
        #self.make_junction(auto_adjust)


    #
    # Aux. functions
    #
    def set_electrode1(self, elec1): 
        if isinstance(elec1, AtomsSystem): self._elec1 = elec1

    def set_electrode2(self, elec2): 
        if isinstance(elec2, AtomsSystem): self._elec2 = elec2

    def set_channel(self, chan):
        if isinstance(chan, AtomsSystem): self._chan = chan

    def get_electrode1(self): return self._electrode1
    def get_electrode2(self): return self._electrode2
    def get_channel(self):   return self._channel


#2345678901234567890123456789012345678901234567890123456789012345678901234567890
    def is_possible(self, tol_angle=10**-6, tol_length=0.1):

        """
Check the given electrodes and channel for possiblity of transport calculation.
        """

        # messages
        info1 = """
Conventionally the transport direction is regarded as the third cell vector and 
perpendicular to the z-direction. 

Currently, the transport direction is fixed to the third cell vector. 
        """
        print (info1)

        # get cell vectors from electrode & channel
        a = self._elec1.get_cell()
        b = self._chan.get_cell()
        ta1 = Vector(a[2])                      # transport dir (elec)
        ta2 = Vector(a[0]); ta3 = Vector(a[1])  # the other dirs (elec)
        tb1 = Vector(b[2])                      # transport dir (chan)
        tb2 = Vector(b[0]); tb3 = Vector(b[1])  # the other dirs (chan)

        # condition 1: trans dir should be orthogonal to the other dirs.
        cond1 = ta1.dot(ta2); cond2 = ta1.dot(ta3)
        cond3 = tb1.dot(tb2); cond4 = tb1.dot(tb3)
        print ("Condition 1: trans dir should be orthogonal to the other directions.")
        print ("dot product between ta1 and ta2 :", cond1)
        print ("dot product between ta1 and ta3 :", cond2)
        print ("dot product between tb1 and tb2 :", cond3)
        print ("dot product between tb1 and tb3 :", cond4)
        print ("\n")

        for cond in [cond1, cond2, cond3, cond4]:
            if (cond > tol_angle):
                print ("--> Cell shape is not orthorhombic.")
                print ("\n")
                return False

        # condition 2: The angle between the other dirs should be same.
        r1 = ta2.angle(ta3)*180./np.pi
        r2 = tb2.angle(ta3)*180./np.pi
        print ("Condition 2: The angle between the other dirs should be same.")
        print ("The angle between ta2 and ta3 =", r1)
        print ("The angle between tb2 and tb3 =", r2)
        print ("\n")

        if r1-r2 > tol_angle:
            print ("--> Cell shape is not matched.")
            print ("\n")
            return False

        # condition 3: The length of the other dirs should be same. 
        l1 = ta2.length(); l2 = ta3.length()
        l3 = tb2.length(); l4 = tb3.length()
        print ("|ta1| and |ta2| :", l1, l2)
        print ("|tb1| and |tb2| :", l3, l4)
        print ("\n")

        cond1 = abs(l1-l3) ; cond2 = abs(l2-l4)
        if (cond1 > tol_length) or (cond2 > tol_length):
            cond1 = abs(l1-l4) ; cond2 = abs(l2-l3)
            if (cond1 > tol_length) or (cond2 > tol_length):
                print ("--> Cell size is not matched.")
                print ("\n")
                return False

        # Final determination
        print ("--> Possible to calculate quantum transport.")
        print ("\n")
        return True


#2345678901234567890123456789012345678901234567890123456789012345678901234567890
    def determine_nkpt(self, elec, max_nk=10, init_nk=1, step_nk=1, opts=None):

        """
obj.determine_nkpt(elec, max_nk=10, init_nk=1, step_nk=1, opts=None)
    
    Automation tool for determination of the number of k-points for electrodes
    (NOT used by users)

    Parameters
    ----------
    elec : AtomsSystem
        an atomic structure as an electrode part

    Optional parameters
    -------------------
    max_nk : integer
        the maximum number of k-points along the shortest lattice vector 
        (longest lattice vector in reciprocal space)
        --> the number of k-points for the other lattice vectors will be 
            automatically determined by their length of reciprocal lattice vector.
    init_nk : integer
        the initial number of k-points along the shortest lattice vector 
        (longest in reciprocal space)
    step_nk : integer
        the increase of k-points 
    opts : dict
        options for SIESTA (refer to siesta.Siesta.get_options().)

    Examples
    --------
    
        """

        from math import ceil

        # cell vector lengths
        cell = elec.get_cell()
        l1 = Vector(cell[0]).length(); l2 = Vector(cell[1]).length(); l3 = Vector(cell[2]).length()
        print ("The length of each cell vectors:")
        print ("%10.6f %10.6f %10.6f\n" % (l1, l2, l3))

        r1 = 1./l1; r2 = 1./l2; r3 = 1./l3
        print ("The 1/length of each cell vectors:")
        print ("%10.6f %10.6f %10.6f\n" % (r1, r2, r3))

        # imposed PBC
        pbc = elec.get_pbc()
        print ("The imposed periodic boundary condition is")
        print (pbc, "along each cell vectors.\n")

        # Check for the periodicity along the transport direction
        if not pbc[2]:
            print ("|========================================|")
            print ("|          CHECK FOR PERIODICITY         |")
            print ("|========================================|")
            print ("Your cell vector is")
            print ("v1 =", "%10.6f %10.6f %10.6f" % tuple(cell[0]))
            print ("v2 =", "%10.6f %10.6f %10.6f" % tuple(cell[1]))
            print ("v3 =", "%10.6f %10.6f %10.6f\n" % tuple(cell[2]))
            print ("At least, the third cell vector (v3) might be periodic.")
            print ("Cannot determine the number of k-points.")
            print ("If your electrode is periodic, please use set_pbc() method.\n")
            raise ValueError("The transport direction (z) should be periodic.")

        # Consider periodic directions only
        rm_ = []; i = 0
        for r in [r1,r2,r3]:
            if elec.get_pbc()[i]: rm_.append(r)
            i += 1
        print (rm_)

        # Find the longest reciprocal lattice vector
        rm = max(rm_)
        r1 = r1/rm; r2 = r2/rm; r3 = r3/rm
        print ("The ratio devided by the maximum 1/length value:")
        print ("%10.6f %10.6f %10.6f\n" % (r1, r2, r3))

        #if not np.array(pbc).any():
        #    print "Your cell vector is"
        #    print "v1 =", "%10.6f %10.6f %10.6f" % tuple(cell[0])
        #    print "v2 =", "%10.6f %10.6f %10.6f" % tuple(cell[1])
        #    print "v3 =", "%10.6f %10.6f %10.6f\n" % tuple(cell[2])
        #    print "At least, the third cell vector (v3) might be periodic."
        #    print "Cannot determine the number of k-points."
        #    print "If your electrode is periodic, please use set_pbc() method."
        #    #print elec.set_pbc.__doc__
        #    #return

        # set kgrid points
        kgrids = []
        i = init_nk
        while i <= max_nk:
            k1, k2, k3 = r1*i, r2*i, r3*i
            k1, k2, k3 = int(ceil(k1)), int(ceil(k2)), int(ceil(k3))
            if not pbc[0]: k1 = 1
            if not pbc[1]: k2 = 1
            if not pbc[2]: k3 = 1
            if [k1,k2,k3] not in kgrids: kgrids.append([k1,k2,k3])
            i += step_nk
        print kgrids
        #import sys
        #sys.exit()

        # run sims over sets of kgrid
        energies = []
        for kgrid in kgrids:

            # simulation object
            sim = s2.Siesta(elec)

            # set simulation options
            sim.set_option('kgrid', kgrid)
            if opts:
                for opt, arg in opts:
                    if opt != 'kgrid': sim.set_option(opt, arg)

            # run siesta
            sim.run()
            e = s2.get_total_energy()
            energies.append(e)

        return kgrids, energies


    def determine_contact_distance(self, elec, chan, opts=None):

        # measure zmin/max of electrode & channel
        elec_ = elec.copy()
        elec_.select_all()
        elec_.sort('z')
        zelec_min = elec_[0][2]
        zelec_max = elec_[-1][2]

        chan_ = chan.copy()
        chan_.select_all()
        chan_.sort('z')
        zchan_min = chan_[0][2]
        zchan_max = chan_[-1][2]

        return



    def make_junction(self, auto_adjust=1): return

    #
    # Main
    #
    def divide_LCR(self): return
