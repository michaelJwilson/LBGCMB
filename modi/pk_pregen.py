#!usr/bin/env python
#
# Some code primarily for playing around for now.
# Makes a class which combines the pre-computed P(k) tables
# to return a single P(k) given the parameters.
#
from __future__ import print_function,division

import numpy             as np
import modi_tools        as tools

from   scipy.interpolate import InterpolatedUnivariateSpline as interpolate


def cleft(lpt,b1,b2,bs,bn,alpha,sn,auto):
    """  
    Returns P(k) given the parameters for the array
    """
    # The column order is:
    # k,P_Z,P_A,P_W,P_d,P_{dd},P_{d^2},
    # P_{d^2d^2},P_{dd^2},P_{s^2},P_{ds^2},
    # P_{d^2s^2},P_{s^2s^2},P_{D2d},P_{dD2d}
    if auto:
        par = np.array([0,1.,1.,1.,b1,b1*b1,b2,b2*b2,b1*b2,\
                            bs,b1*bs,b2*bs,bs*bs,bn,b1*bn])
        par[7] /= 4. #Bad normalization in ZV's files
    else:
        par = np.array([0,1.,1.,1.,b1/2.,0,b2/2.,0,0,bs/2.,0,0,0,bn/2.,0])
    tmp = np.dot(lpt,par)
    tmp-= alpha/2.*lpt[:,0]**2*lpt[:,1]
    return( (lpt[:,0],tmp+sn) )

def cleft_der(pktable,b1,b2,bs,bn,alpha,sn, pnum, auto):
    """  
    Returns the derivative of P(k) from pktable 
    w.r.t. parameter number pnum, with the "usual" order 
    0: b1, 1:b2, 2:bs, 3:bn, 4:alpha, 5:sn, at the 
    position in parameter space passed at "p"

    """
    # The column order is:
    # k,P_Z,P_A,P_W,P_d,P_{dd},P_{d^2},
    # P_{d^2d^2},P_{dd^2},P_{s^2},P_{ds^2},
    # P_{d^2s^2},P_{s^2s^2},P_{D2d},P_{dD2d}
    if auto:
        if pnum==0:     # b1
            dP = pktable[:,4] + 2*b1*pktable[:, 5] \
                 + b2*pktable[:, 8] + bs*pktable[:, 10] + bn*pktable[:, 14]
        elif pnum==1:	# b2
            dP = pktable[:,6] + 2*b2*pktable[:, 7] / 4. 
            + b1*pktable[:,8] + bs*pktable[:, 11]
        elif pnum==2:	# bs
            dP = pktable[:,9] + 2*bs*pktable[:, 12] \
                 + b1*pktable[:, 10]  + b2*pktable[:, 11] 
        elif pnum==3:	# bn
            dP = pktable[:,13] + b1*pktable[:, 14] 
        elif pnum==4:	# alpha
            dP =-0.5*pktable[:,0]**2*pktable[:,1]
        elif pnum==5:	# sn
            dP = 1.0*np.ones_like(pktable[:,0])
        else:
            raise RuntimeError("Unknown pnum="+str(pnum))

    else:
        if pnum==0:		# b1
            dP = 0.5*pktable[:,4]
        elif pnum==1:	# b2
            dP = 0.5*pktable[:,6]
        elif pnum==2:	# bs
            dP = 0.5*pktable[:,9]
        elif pnum==3:	# bn
            dP = 0.5*pktable[:,13]
        elif pnum==4:	# alpha
            dP =-0.5*pktable[:,0]**2*pktable[:,1]
        elif pnum==5:	# sn
            dP = 1.0*np.ones_like(pktable[:,0])
        else:
            raise RuntimeError("Unknown pnum="+str(pnum))
    return( (pktable[:,0],dP) )


class PkCLEFT:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __call__(self,b1,b2,bs,bn,alpha,sn,auto=False):
        """  
        Returns P(k) given the parameters.
        """
        return cleft(self.pktable,b1,b2,bs,bn,alpha,sn,auto)


    def deriv(self, b1,b2,bs,bn,alpha,sn, pnum, auto=False):
        """    
        Returns the derivative of P(k) w.r.t. parameter number pnum, with
        the "usual" order 0: b1, 1:b2, 2:bs, 3:bn, 4:alpha, 5:sn, at the 
        position in parameter space passed at "p"
        """
        return cleft_der(self.pktable, b1, b2, bs, bn, alpha, sn, pnum, auto)
#
    def __init__(self,fn=None, db = None, z = None):
        """    
        Initializes the class by reading in the pre-computed P(k) table.
        The power spectrum table filename is passed as a parameter, if
        nothing is passed a DarkSky cosmology at z=1 is used.
        """
        # Read the data from pre-computed file.
        if fn==None:            
            if db is None:
                db = "../data/RunPB/" # The directory where the data live.
            if z is None:
                z = 2
            self.pktable=np.loadtxt(db +"PTdata/ps00_hh_RunPB_46_z%03d.dat"%(z*100))
        else:
            self.pktable=np.loadtxt(fn)

        #


class PkZel:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __call__(self,b1,b2,sn,auto=False):
        """  
        Returns P(k) given the parameters.
        """
        bs, bn, alpha = 0, 0, 0 #SET THIS
        return cleft(self.pktable,b1,b2,bs,bn,alpha,sn,auto)


#
    def __init__(self,fn=None, db = None, z = None):
        """    
        Initializes the class by reading in the pre-computed P(k) table.
        The power spectrum table filename is passed as a parameter, if
        nothing is passed a DarkSky cosmology at z=1 is used.
        """
        # Read the data from pre-computed file.
        if fn==None:            
            if db is None:
                db = "../data/RunPB/" # The directory where the data live.
            if z is None:
                z = 2
            self.pktable=np.loadtxt(db +"PTdata/ps00_zel_hh_RunPB_46_z%03d.dat"%(z*100))
        else:
            self.pktable=np.loadtxt(fn)

        #

class PkZeft:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __call__(self,b1,b2,alpha,auto=False):
        """  
        Returns P(k) given the parameters.
        """
        bs, bn, sn = 0, 0, 0 #SET THIS
        return cleft(self.pktable,b1,b2,bs,bn,alpha,sn,auto)

#
    def __init__(self,fn=None, db = None, z = None):
        """    
        Initializes the class by reading in the pre-computed P(k) table.
        The power spectrum table filename is passed as a parameter, if
        nothing is passed a DarkSky cosmology at z=1 is used.
        """
        # Read the data from pre-computed file.
        if fn==None:            
            if db is None:
                db = "../data/RunPB/" # The directory where the data live.
            if z is None:
                z = 2
            self.pktable=np.loadtxt(db +"PTdata/ps00_zel_hh_RunPB_46_z%03d.dat"%(z*100))
        else:
            self.pktable=np.loadtxt(fn)

        #


class PkHZeft:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __call__(self,b1,b2,alpha,sn,auto=False):
        """  
        Returns P(k) given the parameters.
        """
        bs, bn = 0,  0 #SET THIS
        return cleft(self.pktable,b1,b2,bs,bn,alpha,sn,auto)

#
    def __init__(self,fn=None, db = None, z = None):
        """    
        Initializes the class by reading in the pre-computed P(k) table.
        The power spectrum table filename is passed as a parameter, if
        nothing is passed a DarkSky cosmology at z=1 is used.
        """
        # Read the data from pre-computed file.
        if fn==None:            
            if db is None:
                db = "../data/RunPB/" # The directory where the data live.
            if z is None:
                z = 2
            self.pktable=np.loadtxt(db +"PTdata/ps00_zel_hh_RunPB_46_z%03d.dat"%(z*100))
        else:
            self.pktable=np.loadtxt(fn)

        #



class Pk8CLEFT(PkCLEFT):
    #### DOES NOT SUPPORT bn ####
    '''Given a pat h (db), redshift(z), varioation(vary) and CLEFT fit params(p)
    interpolate over sigma8 assuming the fodler tree-
    db/ps_response/response_files
    where the template of file names is -
    basefile = "db/ps00_hh_RunPB_46_z%03d.dat"%iz
    responsefile = "db/ps_response/ps00_hh_RunPB_46_z%03d_rp+%02d.dat"%(iz, vary)
    '''
    

    def __init__(self, db =  "../data/RunPB/",  z= 2, vary = [3, 5]):
        '''
        '''
        self.db = db
        self.dbpt = db + 'PTdata/'
        self.z = z
        self.iz = z*100
        self.fn = self.dbpt +"ps00_hh_RunPB_46_z%03d.dat"%self.iz
        self.lin = np.loadtxt(db +'class_output/RunPB00_z%03d_pk.dat'%self.iz)
        
        self.sig80 = 0.8195#tools.sig8(self.lin[:, 0], self.lin[:, 1])
        ###super().__init__(self.fn)
        self.lpt = PkCLEFT(fn = self.fn)
        self.vary = np.array(vary)
        self.vary.sort()
        self.sig8s = self._sig8s()
        self.rsp = self._readfiles()
        self.intpf = self._interp()
        self.pardict = {'s8':-1, 'b1':0, 'b2':1, 'bs2':2, 'bn':3, 'alpha':4, 'sn':5}

    def __call__(self, par ,auto=False):
        """  
        Returns P(k) given the parameters.
        """
        if len(par) == 2:
            p = [0, 0, 0, 0] + list(par) +[auto]
            return self.lpt(*p)
        elif len(par) == 5:
            p = list(par[:3]) +[0] +  list(par[3:])+ [auto]
            return self.lpt(*p)
        elif len(par) == 6:
            s8 = par[0]
            p = list(par[1:4]) + [0] + list(par[4:]) + [auto]
            #rsp = np.zeros_like(self.intpf, dtype = 'f8')
            #for foo in range(rsp.shape[1]):
            #    rsp[:, foo] = np.array([f(s8) for f in self.intpf[:, foo]])
            wts = tools.lagrange2(s8, self.sig8s)
            rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
            return cleft(rsp, *p)

    def deriv(self, par, name, auto=False):
        """    
        Returns the derivative of P(k) w.r.t. parameter number pnum, with
        the "usual" order 0: s8, 1: b1, 2:b2, 3:bs, 4:bn, 5:alpha, 6:sn, at the 
        position in parameter space passed at "p"
        """
        pnum = self.pardict[name]
        if pnum < 0:
            if len(par) == 6:
                s8 = par[0]
                p = list(par[1:4]) + [0] + list(par[4:])
                x1, x2, y = s8*1.01, s8*0.99, []
                for x in [x1, x2]:
                    wts = tools.lagrange2(x, self.sig8s)
                    rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
                    y.append(cleft(rsp, *p, auto = auto))
                der = (y[0][1] - y[1][1])/(x1 - x2)
                return (y[0][0], der)

            else:
                print('For sigma8, parameter length should be 6')

        if pnum >= 0:            
            if len(par) == 2:
                p = [0, 0, 0, 0] + list(par)
                return self.lpt.deriv(*p, pnum = pnum, auto = auto)
            elif len(par) == 5:
                p = list(par[:3]) +[0] +  list(par[3:])
                return self.lpt.deriv(*p, pnum = pnum, auto= auto)
            elif len(par) == 6:
                s8 = par[0]
                p = list(par[1:4]) + [0] + list(par[4:])
                wts = tools.lagrange2(s8, self.sig8s)
                rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
                return cleft_der(rsp, *p, pnum = pnum, auto = auto)



    def _readfiles(self):
        '''Read the response files'''
        files = []
        db = self.dbpt + 'ps_response/'
        for i in self.vary[::-1]:
            files.append(np.loadtxt(db + "ps00_hh_RunPB_46_z%03d_rp-%02d.dat"%(self.iz, i)))
        files.append(np.loadtxt(self.fn))
        for i in self.vary:
            files.append(np.loadtxt(db + "ps00_hh_RunPB_46_z%03d_rp+%02d.dat"%(self.iz, i)))
        return np.array(files)


    def _sig8s(self):
        '''For given 'vary', create array of sigma8s'''
        sig8s = []
        for i in self.vary[::-1]:
            sig8s.append(1 - i/100.)
        sig8s.append(1)
        for i in self.vary:
            sig8s.append(1 + i/100.)
        return np.array(sig8s)*self.sig80

    def _interp(self):
        '''Create interpolation function array to go from sigma8 to reponse array'''
        intpf = []
        rsp = self.rsp
        for foo in range(rsp.shape[1]):
            for boo in range(rsp.shape[2]):
                intpf.append(interpolate(self.sig8s, rsp[:, foo, boo]))
        
        intpf = (np.array(intpf)).reshape(rsp[0].shape)
        return intpf


class Pk8Zel:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __init__(self, db =  "../data/RunPB/",  z= 2, vary = [3, 5]):
        '''
        '''
        self.db = db
        self.dbpt = db + 'PTdata/'
        self.z = z
        self.iz = z*100
        self.fn = self.dbpt +"ps00_zel_hh_RunPB_46_z%03d.dat"%self.iz
        self.lin = np.loadtxt(db +'class_output/RunPB00_z%03d_pk.dat'%self.iz)
        self.signl = (np.trapz(self.lin[:, 1], self.lin[:, 0])/6/np.pi**2)**0.5
        print("0.25*Signl =", 0.25*self.signl)
        self.sig80 = 0.8195
        self.zel = PkZel(fn = self.fn)
        self.vary = np.array(vary)
        self.vary.sort()
        self.sig8s = self._sig8s()
        self.rsp = self._readfiles()
        self.intpf = self._interp()

    def __call__(self, par ,auto=False, pade = False):
        """  
        Returns P(k) given the parameters.
        parameters are s8., b1, b2, sn
        """
        if len(par) == 3:
            if pade:
                sn = par[-1]
                k, p = self.zel(*par, auto = auto)
                p -= sn/(1 + (k*self.signl)**2)
                return (k, p)
            else:
                return self.zel(*par, auto = auto)
        elif len(par) == 4:
            s8 = par[0]
            b1, b2, sn = par[1:]            
            bs, bn, alpha = 0, 0, 0 #SET THIS
            wts = tools.lagrange2(s8, self.sig8s)
            rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
            if pade:
                k, p = cleft(rsp,b1,b2,bs,bn,alpha,sn,auto)
                p -= sn/(1 + (k*self.signl)**2)
                return (k, p)
            else:
                return cleft(rsp,b1,b2,bs,bn,alpha,sn,auto)


    def _readfiles(self):
        '''Read the response files'''
        files = []
        db = self.dbpt + 'ps_response/'
        for i in self.vary[::-1]:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp-%02d.dat"%(self.iz, i)))
        files.append(np.loadtxt(self.fn))
        for i in self.vary:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp+%02d.dat"%(self.iz, i)))
        return np.array(files)

    def _sig8s(self):
        '''For given 'vary', create array of sigma8s'''
        sig8s = []
        for i in self.vary[::-1]:
            sig8s.append(1 - i/100.)
        sig8s.append(1)
        for i in self.vary:
            sig8s.append(1 + i/100.)
        return np.array(sig8s)*self.sig80

    def _interp(self):
        '''Create interpolation function array to go from sigma8 to reponse array'''
        intpf = []
        rsp = self.rsp
        for foo in range(rsp.shape[1]):
            for boo in range(rsp.shape[2]):
                intpf.append(interpolate(self.sig8s, rsp[:, foo, boo]))
        
        intpf = (np.array(intpf)).reshape(rsp[0].shape)
        return intpf




class Pk8Zeft:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __init__(self, db =  "../data/RunPB/",  z= 2, vary = [3, 5]):
        '''
        '''
        self.db = db
        self.dbpt = db + 'PTdata/'
        self.z = z
        self.iz = z*100
        self.fn = self.dbpt +"ps00_zel_hh_RunPB_46_z%03d.dat"%self.iz
        self.lin = np.loadtxt(db +'class_output/RunPB00_z%03d_pk.dat'%self.iz)
        
        self.sig80 = 0.8195
        self.zeft = PkZeft(fn = self.fn)
        self.vary = np.array(vary)
        self.vary.sort()
        self.sig8s = self._sig8s()
        self.rsp = self._readfiles()
        self.intpf = self._interp()

    def __call__(self, par ,auto=False):
        """  
        Returns P(k) given the parameters.
        parameters are s8., b1, b2, sn
        """
        if len(par) == 3:
            return self.zeft(*par, auto = auto)
        elif len(par) == 4:
            s8 = par[0]
            b1, b2, alpha = par[1:]            
            bs, bn, sn = 0, 0, 0 #SET THIS
            wts = tools.lagrange2(s8, self.sig8s)
            rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
            return cleft(rsp,b1,b2,bs,bn,alpha,sn,auto)


    def _readfiles(self):
        '''Read the response files'''
        files = []
        db = self.dbpt + 'ps_response/'
        for i in self.vary[::-1]:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp-%02d.dat"%(self.iz, i)))
        files.append(np.loadtxt(self.fn))
        for i in self.vary:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp+%02d.dat"%(self.iz, i)))
        return np.array(files)

    def _sig8s(self):
        '''For given 'vary', create array of sigma8s'''
        sig8s = []
        for i in self.vary[::-1]:
            sig8s.append(1 - i/100.)
        sig8s.append(1)
        for i in self.vary:
            sig8s.append(1 + i/100.)
        return np.array(sig8s)*self.sig80

    def _interp(self):
        '''Create interpolation function array to go from sigma8 to reponse array'''
        intpf = []
        rsp = self.rsp
        for foo in range(rsp.shape[1]):
            for boo in range(rsp.shape[2]):
                intpf.append(interpolate(self.sig8s, rsp[:, foo, boo]))
        
        intpf = (np.array(intpf)).reshape(rsp[0].shape)
        return intpf



class Pk8HZeft:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __init__(self, db =  "../data/RunPB/",  z= 2, vary = [3, 5]):
        '''
        '''
        self.db = db
        self.dbpt = db + 'PTdata/'
        self.z = z
        self.iz = z*100
        self.fn = self.dbpt +"ps00_zel_hh_RunPB_46_z%03d.dat"%self.iz
        self.lin = np.loadtxt(db +'class_output/RunPB00_z%03d_pk.dat'%self.iz)
        
        self.sig80 = 0.8195
        self.hzeft = PkHZeft(fn = self.fn)
        self.vary = np.array(vary)
        self.vary.sort()
        self.sig8s = self._sig8s()
        self.rsp = self._readfiles()
        self.intpf = self._interp()

    def __call__(self, par ,auto=False):
        """  
        Returns P(k) given the parameters.
        parameters are s8., b1, b2, alpha, sn
        """
        if len(par) == 4:
            return self.hzeft(*par, auto = auto)
        elif len(par) == 5:
            s8 = par[0]
            b1, b2, alpha, sn= par[1:]            
            bs, bn = 0, 0 #SET THIS
            wts = tools.lagrange2(s8, self.sig8s)
            rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
            return cleft(rsp,b1,b2,bs,bn,alpha,sn,auto)


    def _readfiles(self):
        '''Read the response files'''
        files = []
        db = self.dbpt + 'ps_response/'
        for i in self.vary[::-1]:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp-%02d.dat"%(self.iz, i)))
        files.append(np.loadtxt(self.fn))
        for i in self.vary:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp+%02d.dat"%(self.iz, i)))
        return np.array(files)

    def _sig8s(self):
        '''For given 'vary', create array of sigma8s'''
        sig8s = []
        for i in self.vary[::-1]:
            sig8s.append(1 - i/100.)
        sig8s.append(1)
        for i in self.vary:
            sig8s.append(1 + i/100.)
        return np.array(sig8s)*self.sig80

    def _interp(self):
        '''Create interpolation function array to go from sigma8 to reponse array'''
        intpf = []
        rsp = self.rsp
        for foo in range(rsp.shape[1]):
            for boo in range(rsp.shape[2]):
                intpf.append(interpolate(self.sig8s, rsp[:, foo, boo]))
        
        intpf = (np.array(intpf)).reshape(rsp[0].shape)
        return intpf


class Pk8ZA:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __init__(self, db =  "../data/RunPB/",  z= 2, vary = [3, 5]):
        '''
        '''
        self.db = db
        self.dbpt = db + 'PTdata/'
        self.z = z
        self.iz = z*100
        self.fn = self.dbpt +"ps00_zel_hh_RunPB_46_z%03d.dat"%self.iz
        self.lin = np.loadtxt(db +'class_output/RunPB00_z%03d_pk.dat'%self.iz)
        self.signl = (np.trapz(self.lin[:, 1], self.lin[:, 0])/6/np.pi**2)**0.5
        #print("0.25*Signl =", 0.25*self.signl)
        self.sig80 = 0.8195
        self.k = np.loadtxt(self.fn)[:, 0]
        self.pza = np.loadtxt(self.fn)[:, 1]
        self.vary = np.array(vary)
        self.vary.sort()
        self.sig8s = self._sig8s()
        self.rsp = self._readfiles()
        self.intpf = self._interp()

    def __call__(self, par ,auto=False, pade = False):
        """  
        Returns P(k) given the parameters.
        parameters are b1, b11, sn
        """
        if len(par) == 3:
            b1, b11, sn = par[:]
            k, p = self.k, self.pza
        elif len(par) == 4:
            s8, b1, b11, sn = par[:]
            wts = tools.lagrange2(s8, self.sig8s)
            rsp = (wts.reshape(-1, 1, 1)*self.rsp).sum(axis = 0)
            k, p = rsp[:, 0], rsp[:, 1]
        p *= (b1 + b11*k**2)
        if auto:
            p *= (b1 + b11*k**2)
        p += sn
        if pade:
            p -= sn/(1 + (k*self.signl)**2)
            return (k, p)
        else:
            return (k, p)


    def _readfiles(self):
        '''Read the response files'''
        files = []
        db = self.dbpt + 'ps_response/'
        for i in self.vary[::-1]:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp-%02d.dat"%(self.iz, i))[:, :2])
        files.append(np.loadtxt(self.fn)[:, :2])
        for i in self.vary:
            files.append(np.loadtxt(db + "ps00_zel_hh_RunPB_46_z%03d_rp+%02d.dat"%(self.iz, i))[:, :2])
        return np.array(files)

    def _sig8s(self):
        '''For given 'vary', create array of sigma8s'''
        sig8s = []
        for i in self.vary[::-1]:
            sig8s.append(1 - i/100.)
        sig8s.append(1)
        for i in self.vary:
            sig8s.append(1 + i/100.)
        return np.array(sig8s)*self.sig80

    def _interp(self):
        '''Create interpolation function array to go from sigma8 to reponse array'''
        intpf = []
        rsp = self.rsp
        for foo in range(rsp.shape[1]):
            for boo in range(rsp.shape[2]):
                intpf.append(interpolate(self.sig8s, rsp[:, foo, boo]))
        
        intpf = (np.array(intpf)).reshape(rsp[0].shape)
        return intpf


class Pk8Lin:
    """   
    A class which holds the pre-computed P(k) table and has a single
    method which returns the prediction gives the bias parameters.
    """
    def __init__(self, db =  "../data/RunPB/",  z= 2, vary = [3, 5]):
        '''
        '''
        self.db = db
        self.dbpt = db + 'PTdata/'
        self.z = z
        self.iz = z*100
        self.lin = np.loadtxt(db +'class_output/RunPB00_z%03d_pk.dat'%self.iz)
        self.k = self.lin[:, 0]
        self.p = self.lin[:, 1]
        self.signl = (np.trapz(self.lin[:, 1], self.lin[:, 0])/6/np.pi**2)**0.5
        #print("0.25*Signl =", 0.25*self.signl)
        self.sig80 = 0.8195

    def __call__(self, par ,auto=False, pade = False):
        """  
        Returns P(k) given the parameters.
        parameters are b1, b11, sn
        """
        if len(par) == 3:
            b1, b11, sn = par[:]
            k, p = self.k, self.p
        elif len(par) == 4:
            s8, b1, b11, sn = par[:]
            k, p = self.k, self.p*s8**2/self.sig80**2
        p *= (b1 + b11*k**2)
        if auto:
            p *= (b1 + b11*k**2)
        p += sn
        if pade:
            p -= sn/(1 + (k*self.signl)**2)
            return (k, p)
        else:
            return (k, p)


    
class Pk8Hfit():
    #Does not support bn yet!!!!
    '''Lagrange interpolation to get halofit spectra, given a sigma8
    '''
    
    def __init__(self, db =  "../data/RunPB/",  z = 2 ):
        '''
        '''
        self.db = db
        self.z = z
        self.iz = z*100
        self.fn = self.db +"hfit_s8rsp_z%03d.txt"%self.iz
        self.klin = np.loadtxt(self.fn)[:, 0]
        self.hfit = np.loadtxt(self.fn)[:, 1:]
        self.htrue = self.hfit[:, 2]
        self.s80 = 0.8195 #tools.sig8(self.lin[:, 0], self.lin[:, 1])
        self.rat = np.array([0.95, 0.97, 1.0, 1.03, 1.05])
        self.s8s = self.rat*self.s80


    def __call__(self, par, auto = False, pade = False):
        """  
        Returns P(k) given the parameters.
        """
        if type(par) == list:
            if len(par) == 3:
                b1, b11, sn = par[:]
                k, p = self.klin, self.htrue
            elif len(par) == 4:
                s8, b1, b11, sn = par[:]
                k = self.klin
                p = (tools.lagrange2(s8, self.s8s)* self.hfit).sum(axis = -1)
            p *= (b1 + b11*k**2)
            if auto:
                p *= (b1 + b11*k**2)
            p += sn
            if pade:
                p -= sn/(1 + (k*self.signl)**2)
                return (k, p)
            else:
                return (k, p)
        else:
            s8 = par
            return  [self.klin, (tools.lagrange2(s8, self.s8s)* self.hfit).sum(axis = -1)]



class PkNbodyZ:
    """
     A class to hold the P(k) files as a function of z so that we can
     interpolate and get P(k,z) for our halo sample.
    """
    def lagrange_spam(self,z):  # Lagrange interpolating coefficients.
        c = np.ones_like(self.za)
        zzdist = self.za[:, None] - self.za[None, :]
        singular = zzdist == 0
        zzdist[singular] = 1.0
        fac = (z - self.za) / zzdist
        fac[singular] = 1.0
        return(fac.prod(axis=-1))
    def __init__(self,db="../data/RunPB/",auto = False):
        self.db = db
        self.iz=[100,150,200,250,300]
        self.za=np.array(self.iz).astype('f8')/100.0
        self.kk=np.loadtxt(self.db + "hm_z100.pkr")[:,0]
        self.pktable=np.zeros( (len(self.kk),len(self.iz)) )
        for i in range(len(self.iz)):
            if auto:
                self.pktable[:,i] = np.loadtxt(self.db+"hh_z%03d.pkr"%self.iz[i])[:,1]
            else:
                self.pktable[:,i] = np.loadtxt(self.db+"hm_z%03d.pkr"%self.iz[i])[:,1]
            self.pktable[:,i]*= 2.*np.pi**2/self.kk**3
    def __call__(self,z=1.0):
        """z is a scalar.  Returns k,P(k,z)."""
        pofz = np.dot(self.pktable,self.lagrange_spam(z))
        return((self.kk,pofz))
        #
    def interpolate(self,kk,zval):
        """kk is a vector of k values, zval is a scalar.  Returns P(kk,zval)"""
        pofz = np.dot(self.pktable,self.lagrange_spam(zval))
        pk   = np.interp(kk,self.kk,pofz)
        return(pk)
        #




class PkNbodyFineZ:
    """
     A class to hold the P(k) files as a function of z so that we can
     interpolate and get P(k,z) for our halo sample.
    """
    def lagrange_spam(self,z):  # Lagrange interpolating coefficients.
        c = np.ones_like(self.za)
        zzdist = self.za[:, None] - self.za[None, :]
        singular = zzdist == 0
        zzdist[singular] = 1.0
        fac = (z - self.za) / zzdist
        fac[singular] = 1.0
        return(fac.prod(axis=-1))
    def __init__(self,db="../data/RunPB/",auto = False):
        self.db = db
        self.iz=[100,150,200,250,300]
        self.za=np.array(self.iz).astype('f8')/100.0
        self.kk=np.loadtxt(self.db + "hm_fine_z100.txt")[:,0]
        self.pktable=np.zeros( (len(self.kk),len(self.iz)) )
        for i in range(len(self.iz)):
            if auto:
                self.pktable[:,i] = np.loadtxt(self.db+"hh_fine_z%03d.txt"%self.iz[i])[:,1]
            else:
                self.pktable[:,i] = np.loadtxt(self.db+"hm_fine_z%03d.txt"%self.iz[i])[:,1]
            #self.pktable[:,i]*= 2.*np.pi**2/self.kk**3
    def __call__(self,z=1.0):
        """z is a scalar.  Returns k,P(k,z)."""
        pofz = np.dot(self.pktable,self.lagrange_spam(z))
        return((self.kk,pofz))
        #
    def interpolate(self,kk,zval):
        """kk is a vector of k values, zval is a scalar.  Returns P(kk,zval)"""
        pofz = np.dot(self.pktable,self.lagrange_spam(zval))
        pk   = np.interp(kk,self.kk,pofz)
        return(pk)
        #







#class Pk8CLEFT(PkCLEFT):
#    #Does not support bn yet!!!!
#    '''Given a pat h (db), redshift(z), varioation(vary) and CLEFT fit params(p)
#    interpolate over sigma8 assuming the fodler tree-
#    db/ps_response/response_files
#    where the template of file names is -
#    basefile = "db/ps00_hh_RunPB_46_z%03d.dat"%iz
#    responsefile = "db/ps_response/ps00_hh_RunPB_46_z%03d_rp+%02d.dat"%(iz, vary)
#    '''
#    
#
#    def __init__(self, db =  "../data/RunPB/PTdata/",  z= 2, vary = [3, 5]):
#        '''
#        '''
#        self.db = db
#        self.z = z
#        self.iz = z*100
#        self.fn = db+"ps00_hh_RunPB_46_z%03d.dat"%self.iz
#        super().__init__(self.fn)
#        self.vary = np.array(vary)
#        self.vary.sort()
#        self.rsp = self._response()
#        self.n = len(self.rsp)
#
#    def __call__(self, par ,auto=False):
#        """  
#        Returns P(k) given the parameters.
#        """
#        if len(par) == 2:
#            p = [0, 0, 0, 0] + list(par) +[auto]
#            return self.rsp[int((self.n-1)/2)](*p)
#        elif len(par) == 5:
#            p = list(par[:3]) +[0] +  list(par[3:])+ [auto]
#            return self.rsp[int((self.n-1)/2)](*p)
#        elif len(par) == 6:
#            p = list(par[1:4]) +[0] +  list(par[4:])+ [auto]
#            s8 = par[0]
#        else:
#            print('Number of paramters must be 2, 5 or 6. bn is not supported')
#
#        rsp = [self.rsp[i](*p) for i in range(len(self.rsp))]
#        k  = rsp[0][0]
#        l = len(rsp)
#        sig8s = [self.sig8(rsp[i][0], rsp[i][1]) for i in range(l)]
#        isig8k = np.array([interpolate(sig8s, [rsp[foo][1][boo] for foo in range(l)]) for boo in range(k.size)])
#        p8 = [isig8k[foo](s8) for foo in range(k.size)]
#        #p8 = np.array([np.interp(s8, sig8s, [rsp[foo][1][boo] for foo in range(l)]) for boo in range(k.size)])
#        return k, np.array(p8)
#
#    def _response(self):
#        db = self.db + 'ps_response/'
#        iz = self.iz
#        vary = self.vary
#        rsp = []
#        for i in vary[::-1]:
#            rsp.append(PkCLEFT(db + "ps00_hh_RunPB_46_z%03d_rp-%02d.dat"%(iz, i)))
#        rsp.append(PkCLEFT(self.fn))
#        for i in vary:
#            rsp.append(PkCLEFT(db + "ps00_hh_RunPB_46_z%03d_rp+%02d.dat"%(iz, i)))
#        return rsp
#
#    def sig8(self, k, p):
#        kr = k*8
#        wt = 3 * (np.sin(kr)/kr - np.cos(kr))/kr**2
#        if wt is 0:
#            wt = 1
#        sig =  np.trapz(p * wt**2. * k**2, k)/2/np.pi**2
#        return sig**0.5
#    
#    

#Copy of Martin's derivative code
#        if auto:
            #dP = numpy.zeros_like(pktable[:, 0])
            #pass

#        else:
#            if pnum==0:		# b1
#                dP = 0.5*self.pktable[:,4]
#            elif pnum==1:	# b2
#                dP = 0.5*self.pktable[:,6]
#            elif pnum==2:	# bs
#                dP = 0.5*self.pktable[:,9]
#            elif pnum==3:	# bn
#                dP = 0.5*self.pktable[:,13]
#            elif pnum==4:	# alpha
#                dP =-0.5*self.pktable[:,0]**2*self.pktable[:,1]
#            elif pnum==5:	# sn
#                dP = 1.0*np.ones_like(self.pktable[:,0])
#            else:
#                raise RuntimeError("Unknown pnum="+str(pnum))
#        return( (self.pktable[:,0],dP) )

