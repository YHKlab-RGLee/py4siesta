#import xml.etree.ElementTree as ET

#tree = ET.parse('siesta.xml')
#root = tree.getroot()


class XmlObject(object):

    def __init__(self, tag='', atr='', txt='', lv=1):
        self._children = []
        self.add_tag(tag)
        self.add_atr(atr)
        self.add_txt(txt)
        self.set_level(lv)
    
    def add_tag(self, tag):
        self._tag = tag

    def add_atr(self, atr):
        self._atr = atr

    def add_txt(self, txt):
        self._txt = txt

    def add_child(self, child):
        self._children.append(child)

    def set_level(self, lv):
        self._level = lv

    def __getitem__(self, i):
        return self._children[i]

    def __len__(self):
        return len(self._children)


class SiestaXmlObject(XmlObject):

    def is_root(self):
        if self._level == 1: return True
        else: return False


    def get_initial_structure(self):

        if self.is_root():

            atoms = []
            for atomobj in self._children[2][0][0]:
                x = float(atomobj._atr['x3']) 
                y = float(atomobj._atr['y3'])
                z = float(atomobj._atr['z3'])
                symb = atomobj._atr['elementType']
                atoms.append( Atom(symb, [x,y,z]) )

            cell = []
            for cellv in self._children[2][1]:
                v1, v2, v3 = cellv._txt.split()
                v1 = float(v1); v2 = float(v2); v3 = float(v3)
                cell.append([v1, v2, v3])

            return AtomsSystem(atoms, cell=np.array(cell)/units.ang2bohr)


    def get_title(self):
        if self.is_root(): return self._children[3][0][0]._txt


    def get_nkpoints(self):
        if self.is_root(): return int(self._children[4][0][0]._txt)

    def get_kpoints(self):
        if self.is_root():
            kpts = []
            kwts = []
            for kpt in  self._children[4][1:]:
                if kpt._tag == "kpoint":
                    kx, ky, kz = kpt._atr['coords'].split()
                    kx = float(kx); ky = float(ky); kz = float(kz)
                    kpts.append([kx,ky,kz])
                    kw = float(kpt._atr['weight'])
                    kwts.append(kw)
            return np.array(kpts), np.array(kwts)


    def get_kdispl(self):
        if self.is_root():
            dkx, dky, dkz = self._children[6][0]._txt.split()
            dkx = float(dkx); dky = float(dky); dkz = float(dkz)
            return np.array([dkx, dky, dkz])


    def get_mdsteps(self):
        if self.is_root():
            MD_steps = {}
            MD_count = 0
            # For all chilren,
            for child in self._children:

                # Find MD module
                if (child._tag == "module") and ('dictRef' in child._atr.keys()):
                    if child._atr['dictRef'] == "MD":
                        MD_count += 1

                        # Inside a MD step...
                        atoms_sw = 0
                        atoms_init = []
                        atoms_fin = []
                        cell_sw = 0
                        cell_init = []
                        cell_fin = []
                        SCF = []
                        E_KS = 0.
                        forces = []

                        for grandchild in child._children:
                            #print grandchild._tag, grandchild._atr, grandchild._txt

                            if grandchild._tag == "molecule":

                                atoms_ = []
                                for atomobj in grandchild[0]:
                                    x = float(atomobj._atr['x3']) 
                                    y = float(atomobj._atr['y3'])
                                    z = float(atomobj._atr['z3'])
                                    symb = atomobj._atr['elementType']
                                    #print atoms_sw, symb, x, y, z
                                    atoms_.append( Atom(symb, [x,y,z]) )

                                atoms_fin = atoms_

                                #if atoms_sw == 0:
                                #    atoms_init = atoms_
                                #    atoms_sw = 1
                                #elif atoms_sw == 1:
                                #    atoms_final = atoms_
                                #    atoms_sw = 2
                                #else: pass

                            elif grandchild._tag == "lattice":

                                cell_ = []
                                for cellv in grandchild:
                                    v1, v2, v3 = cellv._txt.split()
                                    v1 = float(v1); v2 = float(v2); v3 = float(v3)
                                    #print cell_sw, [v1, v2, v3]
                                    cell_.append([v1, v2, v3])

                                cell_fin = cell_

                                #if cell_sw == 0:
                                #    cell_init = cell_
                                #    cell_sw = 1
                                #elif cell_sw == 1:
                                #    cell_final = cell_
                                #    cell_sw = 2
                                #else: pass

                            elif (grandchild._tag == "module") and ('dictRef' in grandchild._atr.keys()):

                                if grandchild._atr['dictRef'] == "SCF" and grandchild._atr['serial'] != "1":
                                    serial_scf = int(grandchild._atr['serial'])
                                    Eharrs = float(grandchild[0][0][0]._txt)
                                    FreeE  = float(grandchild[0][1][0]._txt)
                                    Ef     = float(grandchild[0][2][0]._txt)
                                    #print "SCF", serial_scf, "Eharrs =", Eharrs, "FreeE =", FreeE, "Ef =", Ef
                                    SCF.append( [serial_scf, Eharrs, Ef] )

                            elif (grandchild._tag == "module") and ('title' in grandchild._atr.keys()):

                                if grandchild._atr['title'] == "SCF Finalization":
                                    E_KS     = float(grandchild[0][0][0]._txt)
                                    #E_KS_egg = float(grandchild[0][1][0]._txt)
                                    forces_  = grandchild[1][0][0]._txt.split()
                                    rows = int(grandchild[1][0][0]._atr['rows'])
                                    cols = int(grandchild[1][0][0]._atr['columns'])
                                    forces__ = []
                                    for f in forces_: forces__.append(float(f))
                                    forces__ = np.array(forces__)
                                    forces = forces__.reshape((cols, rows))
                                    #print "forces ="
                                    #print forces

                            else: pass

                        #for at in atoms_init: print at
                        #print cell_init
                        #print AtomsSystem(atoms_init, cell=cell_init)
                        #print SCF
                        #print forces
                        #for at in atoms_fin: print at
                        #print cell_fin

                        #atoms_1 = AtomsSystem(atoms_init, cell=cell_init)
                        atoms_2 = AtomsSystem(atoms_fin, cell=cell_fin)

                        MD_steps[MD_count] = [SCF, E_KS, forces, atoms_2]

            return MD_steps

"""
    def get_final(self):
        if not self.is_root(): return

        # For all chilren,
        for child in self._children:
            if (child._tag == 'module') and ('title' in child._atr.keys()):
                # Finalization part
                if child._atr['title'] == 'SCF Finalization':
                    for grandchild in child:

                        if (grandchild._tag == "propertyList") and ('title' in grandchild._atr.keys())
"""                   


def read_siesta_xml(obj, lv, ith):

    # info
    tag = obj.tag.split('}')[-1]
    atr = obj.attrib
    txt = obj.text
    xmlobj = SiestaXmlObject(tag, atr, txt, lv)

    print "    "*(lv-1), "depth lv =", lv, ith, "-th", \
          "tag =", obj.tag.split('}')[-1], "attrib =", obj.attrib, "text =", obj.text, '\n'

    # children
    i = 1
    for child in obj:
        info1 = read_siesta_xml(child, lv+1, i)
        xmlobj.add_child(info1)
        i += 1

    return xmlobj


#xmlobj = read_siesta_xml(root, 1, 1)
#print xmlobj.get_initial_structure()
#print xmlobj.get_title()
#print xmlobj.get_nkpoints()
#print xmlobj.get_kpoints()
#print xmlobj.get_kdispl()
#dic_md = xmlobj.get_mdsteps()
