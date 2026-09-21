def rotate(self, angle, rot_center, rot_vector):
    """
    rotate a vector using Rodrigues' rotation formula
    >>> angle = 30.0 # in degree
    >>> rot_center = Vector(0,0,0) # center of rotation
    >>> rot_vector = Vector(0,0,1) # rotation axis
    >>> vector.rotate(angle, rot_center, rot_vector)
    """

    rot_v = rot_vector/rot_vector.length()       # rotation axis vector
    obj_v = Vector(self.x(), self.y(), self.z()) # original vector
    ux, uy, uz = rot_v.x(), rot_v.y(), rot_v.z() # rotation axis vector component
    ax, ay, az = obj_v - rot_center              # original vector component centered at (0,0,0)
    ang_rad = angle*pi/180.                      # rotation angle in radian
    cos_a = cos(ang_rad); sin_a = sin(ang_rad)   # cos, sin value of rotation angle
    umat = np.matrix([ax, ay, az]).T             # matrix form of original vector
    # rotation matrix
    Rmat = np.matrix( [[ cos_a + (ux**2)*(1-cos_a),  ux*uy*(1-cos_a) - uz*sin_a, ux*uz*(1-cos_a) + uy*sin_a ], 
                       [ uy*ux*(1-cos_a) + uz*sin_a, cos_a + (uy**2)*(1-cos_a),  uy*uz*(1-cos_a) - ux*sin_a ],
                       [ uz*ux*(1-cos_a) - uy*sin_a, uz*uy*(1-cos_a) + ux*sin_a, cos_a + (uz**2)*(1-cos_a)  ]] )
    # rotated vector componets
    bx, by, bz = Rmat*umat
    # rotated vector vector component centered at the original center
    res_v = Vector(bx, by, bz) + rot_center
    return res_v
