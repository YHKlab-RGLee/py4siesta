from atoms import *

class Transport(object):
    """
Transport(electrode, channel)
    
    Class for representing a quantum transport simulation

    Parameters
    ----------
    electrode : AtomsSystem
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

    __slots__ = ['_electrode1', '_electrode2', '_channel', '_trans_dir_elec', '_trans_dir_chan']

    #
    # Initialization
    #
    def __init__(self, elec1, chan, elec2, trans_dir_elec=3, trans_dir_chan=3, auto_adjust=0):

        # 1. set electrode and channel
        self.set_electrode1(elec1)
        self.set_electrode2(elec2)
        self.set_channel(chan)
        #self.set_trans_dir_electrode(trans_dir_elec)
        #self.set_trans_dir_channel(trans_dir_chan)

        # 2. check: cell shape & size
        self.is_possible(trans_dir_elec, trans_dir_chan)

        # 3. adjust transport direction
        #self.adjust_junction(trans_dir_elec, trans_dir_chan)

        # 4. make junction model
        self.make_junction(auto_adjust)


    #
    # Aux. functions
    #
    def set_electrode1(self, elec1): 
        if isinstance(elec1, AtomsSystem): self._electrode1 = elec1

    def set_electrode2(self, elec2): 
        if isinstance(elec2, AtomsSystem): self._electrode2 = elec2

    def set_channel(self, chan):
        if isinstance(chan, AtomsSystem): self._channel = chan

    def set_trans_dir_electrode(self, trans_dir):
        if trans_dir in [1, 2, 3]: self._trans_dir_elec = trans_dir

    def set_trans_dir_channel(self, trans_dir):
        if trans_dir in [1, 2, 3]: self._trans_dir_elec = trans_dir

    def get_electrode1(self): return self._electrode1
    def get_electrode2(self): return self._electrode2
    def get_channel(self):   return self._channel


    #def is_periodic(self):
    #    p1 = self._electrode.get_pbc()
    #    p2 = self._channel.get_pbc()
    #    d1 = p1[0] and p2[0]
    #    d2 = p1[1] and p2[1]
    #    if d1 or d2: return True
    #    else: return False


    def is_possible(self, trans_dir1=3, trans_dir2=3, tol_angle=10**-6, tol_length=0.1):

        # messages
        warn1 = """
        1. Conventionally the transport direction is regarded as the third cell vector along z-direction. 
        """

        # get cell vectors from electrode & channel
        a = self._electrode1.get_cell()
        b = self._channel.get_cell()
        ia = [1, 2, 3]; ib = [1, 2, 3]
        iia = ia.pop(trans_dir1-1)
        iib = ib.pop(trans_dir2-1)
        ta1 = Vector(a[trans_dir1-1]) # transport dir (elec)
        tb1 = Vector(b[trans_dir2-1]) # transport dir (chan)
        ta2 = Vector(a[ia[0]-1]); ta3 = Vector(a[ia[-1]-1])  # the other dirs (elec)
        tb2 = Vector(b[ib[0]-1]); tb3 = Vector(b[ib[-1]-1])  # the other dirs (chan)
        #print "Transport dir for electrode is %i:" % iia, "ta1 =", ta1
        #print "Transport dir for channel   is %i:" % iib, "tb1 =", tb1
        #print "\n"
        #print "The other vectors for electrode are\n", "ta2 =", ta2, "ta3 =", ta3
        #print "The other vectors for channel   are\n", "tb2 =", tb2, "tb3 =", tb3
        #print "\n"

        # condition 1: trans dir should be orthogonal to the other dirs.
        cond1 = ta1.dot(ta2); cond2 = ta1.dot(ta3)
        cond3 = tb1.dot(tb2); cond4 = tb1.dot(tb3)
        print "dot product between ta1 and ta2 :", cond1
        print "dot product between ta1 and ta3 :", cond2
        print "dot product between tb1 and tb2 :", cond3
        print "dot product between tb1 and tb3 :", cond4
        print "\n"

        for cond in [cond1, cond2, cond3, cond4]:
            if (cond > tol_angle):
                print "--> Cell shape is not orthorhombic."
                print "\n"
                return False

        # condition 2: The angle between the other dirs should be same.
        r1 = ta2.angle(ta3)*180./np.pi
        r2 = tb2.angle(ta3)*180./np.pi
        print "The angle between ta2 and ta3 =", r1
        print "The angle between tb2 and tb3 =", r2
        print "\n"

        if r1-r2 > tol_angle:
            print "--> Cell shape is not matched."
            print "\n"
            return False

        # condition 3: The length of the other dirs should be same. 
        l1 = ta2.length(); l2 = ta3.length()
        l3 = tb2.length(); l4 = tb3.length()
        print "|ta1| and |ta2| :", l1, l2
        print "|tb1| and |tb2| :", l3, l4
        print "\n"

        cond1 = abs(l1-l3) ; cond2 = abs(l2-l4)
        if (cond1 > tol_length) or (cond2 > tol_length):
            cond1 = abs(l1-l4) ; cond2 = abs(l2-l3)
            if (cond1 > tol_length) or (cond2 > tol_length):
                print "--> Cell size is not matched."
                print "\n"
                return False

        # Final determination
        print "--> Possible to calculate quantum transport."
        print "\n"
        return True


    def adjust_junction(self, trans_dir1, trans_dir2, tol=10**-6):

        # 1. define electrode, channel, and z axis
        v1 = Vector(self._electrode.get_cell()[trans_dir1-1])
        v2 = Vector(self._channel.get_cell()[trans_dir2-1])
        v3 = Vector(0,0,1)
        print "Transport dir for electrode is %i:" % trans_dir1, v1
        print "Transport dir for channel   is %i:" % trans_dir2, v2
        print "\n"

        # 2. rotate electrode unitcell --> transport direction is z axis
        elec_ = self._electrode.copy()
        print "Before rotation:"
        print elec_.get_cell()
        rot_angle = np.arccos(v1.dot(v3)/(abs(v1)*abs(v3)))*180./np.pi
        rot_axis = v1.cross(v3)

        if abs(rot_axis) > tol:
            print "--> Adjustment\n"
            elec_.select_all()
            elec_.rotate(rot_angle, Vector(0,0,0), rot_axis, with_cell=1)
            print "rot_angle1 =", rot_angle
            print "rot_axis1  =", rot_axis
            print "After rotation for z:"
            print elec_.get_cell()

            # regrrange
            index = [0, 1, 2]
            cell = elec_.get_cell()
            vt3 = Vector(cell[trans_dir1-1]); index.pop(trans_dir1-1)
            vt2 = Vector(cell[index[-1]])   ; index.pop(-1)
            vt1 = Vector(cell[index[0]])
            cell = [vt1, vt2, vt3]
            elec_.set_cell(cell)
            print "After rearrange:"
            print elec_.get_cell()
      
            # vt1 // x
            rot_angle = np.arccos(vt1.dot(Vector(1,0,0))/(abs(vt1)*abs(Vector(1,0,0))))*180./np.pi       
            rot_axis = Vector(0,0,1)
            print "rot_angle1 =", rot_angle
            print "rot_axis1  =", rot_axis
            elec_.select_all()
            elec_.rotate(rot_angle, Vector(0,0,0), rot_axis, with_cell=1)
            elec_.rotate(180., Vector(0,0,0), Vector(0,0,1), with_cell=1)
            print "After rotation for x:"
            print elec_.get_cell()
            print "\n"
        self.set_electrode(elec_)
      
        # 3. rotate channel unitcell --> transport direction is z axis
        chan_ = self._channel.copy()
        print "Before rotation:"
        print chan_.get_cell()
        rot_angle = np.arccos(v2.dot(v3)/(abs(v2)*abs(v3)))*180./np.pi
        rot_axis = v2.cross(v3)

        if abs(rot_axis) > tol:
            print "--> Adjustment"
            chan_.select_all()
            chan_.rotate(rot_angle, Vector(0,0,0), rot_axis, with_cell=1)
            print "rot_angle1 =", rot_angle
            print "rot_axis1  =", rot_axis
            print "After rotation for z:"
            print chan_.get_cell()

            # regrrange
            index = [0, 1, 2]
            cell = chan_.get_cell()
            vt3 = Vector(cell[trans_dir2-1]); index.pop(trans_dir2-1)
            vt2 = Vector(cell[index[-1]])   ; index.pop(-1)
            vt1 = Vector(cell[index[0]])
            cell = [vt1, vt2, vt3]
            chan_.set_cell(cell)
            print "After rearrange:"
            print chan_.get_cell()
            print "\n"
     
            # vt1 // x
            rot_angle = np.arccos(vt1.dot(Vector(1,0,0))/(abs(vt1)*abs(Vector(1,0,0))))*180./np.pi       
            rot_axis = Vector(0,0,1)
            print "rot_angle1 =", rot_angle
            print "rot_axis1  =", rot_axis
            chan_.select_all()
            chan_.rotate(rot_angle, Vector(0,0,0), rot_axis, with_cell=1)
            chan_.rotate(180., Vector(0,0,0), Vector(0,0,1), with_cell=1)
            print "After rotation for x:"
            print chan_.get_cell()
            print "\n"
        self.set_channel(chan_)
     
        return 
     

    def make_junction(self, auto_adjust=1): return


    def determine_kpts(self): return
    def determine_contact_distance(self): return

    #
    # Main
    #
    def divide_LCR(self): return
