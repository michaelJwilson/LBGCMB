import  glob

import  pylab  as  pl
import  numpy  as  np


files = glob.glob("./*.sed")
files = sorted(files)

for i, file in enumerate(files):    
    ## if file.startswith("../sn-"):
    ##    pass

    data = np.loadtxt(file, comments='#')

    pl.plot(data[:,0], data[:,1])
        
    pl.title(file)

    pl.xlim(4e3, 1.5e4)

    pl.xlabel(r"$\AA$")
    
    pl.show()
        
