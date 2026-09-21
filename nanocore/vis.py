import os, glob, re
from . io import write_xsf, write_pdb




def show_xcrysden(atoms, reset=False, **options):
    
    n = 1
        
    # If "reset", remove all the previous temporary files.
    if reset != False:
        files = glob.glob("temp_xsf/*")
        for i in files:
            os.remove("%s" %i)

    # If the directory "temp" does not exist, make it.
    if os.path.isdir("temp_xsf") == False:
        os.mkdir("temp_xsf")
        
    # Find existing temporary xsf files.
    files = glob.glob("temp_xsf/temp_*.xsf")
    
    # Depending on the exsiting files, determine the file name...
    if files == []:
        pass
    else:
        n2 = []
        for i in files:
            n1 = re.findall('\d+', i)
            n2 += n1
        n2.sort(reverse=True, key=int)
        n= int(n2[0])
            
                
    # Make a xsf file. 
    name = 'temp_xsf/temp_%d.xsf' % n
    write_xsf(name, atoms)

    # Make a xcrysden sript file based on input options.
    write_xsf_script(name, "xsf_script.tcl", **options)

    # Load XcrySDen...
    os.system("xcrysden --script xsf_script.tcl &")
    #os.system("rm -rf temp_xsf")
    #os.system("rm xsf_script.tcl")


def write_xsf_script(xsf_name, tcl_name=None, **options):

    # First, determine script file name...
    default = 'xsf_script.tcl'
    if tcl_name is None: tcl_name = default

    # Write the script file ...
    xs_load(xsf_name,tcl_name)

    # Update the script options...
    if "rot" in options:
        xs_rotate(tcl_name, **options["rot"])
    else:
        xs_rotate(tcl_name, [("x",90),("y",-120),("x",-15)] )

    if "zoom" in options:
        xs_zoom(tcl_name, options["zoom"])
    else:
        xs_zoom(tcl_name, 0.7)
        
    if "full" in options:
        xs_full_screen(tcl_name)

def xs_load(xsf_name,tcl_name):
    f = open(tcl_name, "w")
    f.write("scripting::open --xsf %s\n" % xsf_name)
    f.write("scripting::display on crystal-cells\n")
    f.write("scripting::display on coordinate-system\n")
    f.write("scripting::display on perspective\n")
    f.write("scripting::display as crystal-cells rods\n")
    f.write("scripting::display as cell-unit cell\n")
    #f.write("set myParam(atmCol(1)) {1.0 1.0 1.0}\n")    
    #f.write("scripting::load_myParam\n")
    f.close()

def xs_rotate(tcl_name, rot_opts):
    f = open(tcl_name, 'a')
    print (rot_opts)
    for opt in rot_opts:
        print (opt[0],opt[1])
        f.write("scripting::rotate %s %s\n" % (opt[0],opt[1]))
    f.close()

def xs_zoom(tcl_name, magn):
    f = open(tcl_name, 'a')
    f.write("scripting::zoom %s \n" % magn)
    f.close()

def xs_full_screen(tcl_name):
    f = open(tcl_name, 'a')
    f.write("scripting::displayWindow fullscreen\n")
    f.close()
