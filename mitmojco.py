#!/usr/bin/env python3
#==============================================
#---------- by Dmitry R. Gulevich -------------
#---------- drgulevich@gmail.com --------------
#==============================================
#-----------------------------------------------------------------------------------------
# This is Python mitmojco module for processing tunnel current amplitudes.
# For detailed examples of using the module see Jupyter Python notebook ``amplitudes.ipynb``
#
# tca_bcs(T, Delta1, Delta2): 
#       Calculate of tunnel current amplitudes from BCS expressions.
#
# tca_smbcs(T, Delta1, Delta2, dsm):
#       Smooth Riedel peaks in BCS tunnel current amplitudes. 
#
# new_fit(x, Jpair_data, Jqp_data, maxNterms, thr):
#       Fit tunnel current amplitudes by a sum of complex exponentials.
#
# parsave(pAB,filename):
#       Save fit parameters p,A,B to file.
#
# tca_fit(pAB):
#       Return fitted TCAs for given parameters 
#
# load_fit_parameters(filename):
#       Load TCA fit parameters from file 
#
# load_fit(x, filename):
#       Load TCA fit parameters from file and calculate data points at x
#-----------------------------------------------------------------------------------------
import numpy as np
from scipy.integrate import quad
from scipy.optimize import leastsq
from functools import partial
import time

e0 = 1.6e-19 # electron charge (SI units)
kB = 1.38e-23 # Boltzmann constant (SI units)


#-----------------------------------------------------------------------------------------------
# Function description:
#     Calculation of tunnel current amplitudes 
#     from BCS expressions of Larkin and Ovchinnikov [1]
#     as given by formulas (A5)-(A8) in Appendix of Ref.[2].
#     Tunnel current amplitudes are normalized to Vg/R.
# Input:
#     T -- temperature in K
#     Delta1 -- superconducting gap in meV
#     Delta2 -- superconducting gap in meV
# Output:
#     Repair(x) -- real part of the pair current
#     Impair(x) -- imaginary part of the pair current
#     Reqp(x) -- real part of the quasiparticle current
#     Imqp(x) -- imaginary part of the quasiparticle current
# References:
#     [1] A. I. Larkin and Yu. N. Ovchinnikov, Sov. Phys. JETP 24, 1035 (1967).
#     [2] D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, Phys. Rev. B 96, 024515 (2017).
#-----------------------------------------------------------------------------------------------
def tca_bcs(T, Delta1, Delta2):

    ### define normalized gaps d1 and d2 (0 < d1 <= d2)
    d1=min(Delta1,Delta2)/(Delta1+Delta2)
    d2=1.-d1
    Vg = 1.e-3*(Delta1+Delta2) # gap voltage in Volts
    b=e0*Vg/(2*kB*T) # parameter alpha in Ref. PRB 96, 024515 (2017)

    dd1 = d1*d1
    dd2 = d2*d2
    d21 = d2-d1
    
    ###---------------------------------------------------
    ###          Define tunnel current amplitudes
    ###---------------------------------------------------
    ### In the calculations we assume x>0 without loosing 
    ### the generality (x<0 are obtained by symmetry)

    #------- Re pair current -------
    
    def Repair_integrand_part1(y, x):	
        return np.tanh(b*abs(y))/( np.sqrt(dd1-(y-x)*(y-x)) * np.sqrt(y*y-dd2) )
    
    def Repair_integrand_part2(y, x):
        return np.tanh(b*abs(y))/( np.sqrt(y*y-dd1) * np.sqrt(dd2-(y+x)*(y+x)) )
    
    @np.vectorize
    def Repair(x):
        if(x+d1>d2):
            part1 = quad(Repair_integrand_part1, max(x-d1,d2), x+d1, args=x, limit=200)[0]
        else:
            part1=0.
        if(-x-d2<-d1):    
            part2 = quad(Repair_integrand_part2, -x-d2, min(-x+d2,-d1), args=x, limit=200)[0]
        else:
            part2=0.
        if(-x+d2>d1):
            part3 = quad(Repair_integrand_part2, d1, -x+d2, args=x, limit=200)[0]    
        else:
            part3=0.
        return 0.5*d1*d2*(part1 + part2 + part3)

    #------- Im pair current -------

    def Impair_integrand(y, x):
        return ( np.tanh(b*(y+x)) - np.tanh(b*y) )*np.sign(y)*np.sign(y+x)/( np.sqrt(y*y-dd1) * np.sqrt((y+x)*(y+x)-dd2) )

    @np.vectorize
    def Impair(x):
        part1 = quad(Impair_integrand, -np.inf, min(-d1,-x-d2), args=x, limit=200)[0]
        part2 = quad(Impair_integrand, max(d1,-x+d2), np.inf, args=x, limit=200)[0]
        if(-x+d2<-d1):
            part3 = quad(Impair_integrand, -x+d2, -d1, args=x, limit=200)[0]    
        else:
            part3=0.
        return 0.5*d1*d2*(part1 + part2 + part3)   
    
    #------- Re qp current -------

    def Reqp_integrand_part1(y, x):
        return abs(y)*np.tanh(b*y)*(y-x)/( np.sqrt(y*y-dd1) * np.sqrt(dd2-(y-x)*(y-x)) )

    def Reqp_integrand_part1b_d1_xd2(y, x):
        return abs(y)*np.tanh(b*y)*(y-x)/np.sqrt((y+d1)*(d2+y-x))

    def Reqp_integrand_part1b_xd2_xd2(y, x):
        return abs(y)*np.tanh(b*y)*(y-x)/np.sqrt(y*y-dd1)

    def Reqp_integrand_part2_xd1_d2(y, x):
        return abs(y)*np.tanh(b*y)*(y+x)/np.sqrt((d1-y-x)*(d2-y))

    def Reqp_integrand_part2_xd1_xd1(y, x):
        return abs(y)*np.tanh(b*y)*(y+x)/np.sqrt(y*y-dd2)

    @np.vectorize
    def Reqp(x):
        if(x-d2<-d1):
            part1a = quad(Reqp_integrand_part1, x-d2, -d1, args=x, limit=200)[0]
        else:
            part1a = 0.           
            
        if(x-d2<d1):    
            part1b = quad(Reqp_integrand_part1b_d1_xd2, d1, x+d2, args=x, limit=200, weight='alg', wvar=(-0.5,-0.5))[0]
        else:
            part1b = quad(Reqp_integrand_part1b_xd2_xd2, x-d2, x+d2, args=x, limit=200, weight='alg', wvar=(-0.5,-0.5))[0]
    
        if(-x-d1<-d2):
            if(-d2<-x+d1):
                part2 = quad(Reqp_integrand_part2_xd1_d2, -x-d1, -d2, args=x, limit=200, weight='alg', wvar=(-0.5,-0.5))[0]
            else:
                part2 = quad(Reqp_integrand_part2_xd1_xd1, -x-d1, -x+d1, args=x, limit=200, weight='alg', wvar=(-0.5,-0.5))[0]
        else:
            part2 = 0.
        return -0.5*(part1a + part1b + part2)
    
    #------- Im qp current -------
   
    def Imqp_integrand(y, x):
        return ( np.tanh(b*(y+x))-np.tanh(b*y) )*abs(y)*abs(y+x)/( np.sqrt((y+x)*(y+x)-dd1) * np.sqrt(y*y-dd2) )

    @np.vectorize
    def Imqp(x):
        part1 = quad(Imqp_integrand, -np.inf, min(-d2,-x-d1), args=x, limit=200)[0]
        if(-x+d1<-d2):
            part2 = quad(Imqp_integrand, -x+d1, -d2, args=x, limit=200)[0]
        else:
            part2 = 0.
        part3 = quad(Imqp_integrand, d2, np.inf, args=x, limit=200)[0]
        return 0.5*(part1 + part2 + part3)

    #------- Functions to be returned at arbitrary x (values x<0 obtained by symmetry)

    def Jpair(x):
        absx = np.maximum(np.abs(x),1.e-5)
        return Repair(absx) + np.sign(x)*1j*Impair(absx)

    def Jqp(x):
        absx = np.maximum(np.abs(x),1.e-5)
        return Reqp(absx) + np.sign(x)*1j*Imqp(absx)

    return (Jpair,Jqp)


#-------------------------------------------------------------------------
# Function description: 
#     Smoothing Riedel peaks of bare BCS tunnel current amplitudes (TCAs). 
#     We apply the smoothing procedure as described in Ref. [1-4]. Smoothing parameter 'dsm' 
#     used here is equal to the parameter $\delta$ in [1-4]. Note, that the gaps d1 and d2 
#     are interchanged as compared to Ref.[1]: in our calculations we use 0<d1<d2 as in Ref.[2-4].
# Output:
#     Jpair_smooth(x)
#     Jqp_smooth(x)
# References:
#    [1] A. B. Zorin, I. O. Kulik, K. K. Likharev, and J. R. Schrieffer, 
#        Sov. J. Low Temp. Phys. 5, 537 (1979).
#    [2] D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, Phys. Rev. B 96, 024515 (2017).
#    [3] D. R. Gulevich, V. P. Koshelets, F. V. Kusmartsev, arXiv:1709.04052 (2017).
#    [4] D. R. Gulevich, L. V. Filippenko, V. P. Koshelets, arXiv:1809.01642 (2018).
#-------------------------------------------------------------------------
def tca_smbcs(T, Delta1, Delta2, dsm):

    ### define normalized gaps d1 and d2 (0 < d1 <= d2)
    d1=min(Delta1,Delta2)/(Delta1+Delta2)
    d2=1.-d1
    Vg = 1.e-3*(Delta1+Delta2) # gap voltage in Volts
    b=e0*Vg/(2*kB*T) # parameter alpha in Ref. PRB 96, 024515 (2017)

    ### Supplementary variables
    dd1 = d1*d1
    dd2 = d2*d2
    d21 = d2-d1
    expb = np.exp(b)
    invexpb = 1./expb

#    Repair,Impair,Reqp,Imqp = tca_bcs(T,Delta1,Delta2)
    Jpair,Jqp = tca_bcs(T,Delta1,Delta2)

    IP0=Jpair(0.).real

    ### The smoothing procedure is different for symmetric junction (d1=d2) due 
    ### to the difference d2-d1 appearing in the denominator: the limit d2-d1->0 is 
    ### taken analytically to avoid numerical errors. 
    if(d21<0.001):
        symmetric_junction=True
        print('# Symmetric junction assumed (Delta1=Delta2)')
    else:
        symmetric_junction=False
        print('# Asymmetric junction (Delta1!=Delta2)')

    ### Smoothing for Repair, Reqp
    def dRe(x):
        sqpos = (x-1.)*(x-1.)
        sqneg = (x+1.)*(x+1.)    
        return -x*(1./np.pi)*IP0*0.5*np.log( ((sqpos+dsm*dsm)/sqpos) * (sqneg/(sqneg+dsm*dsm)) )

    ### Smoothing for Repair, -Reqp at x=d2-d1 
    def dRe_minus(x):
        return np.pi*x*np.sqrt(d1*d2)*(np.tanh(b*d2)-np.tanh(b*d1)) * 0.5*( (2./np.pi)*np.arctan((x-d21)/dsm) 
                - np.sign(x-d21) + (2./np.pi)*np.arctan((x+d21)/dsm) - np.sign(x+d21) ) /(4.*d21)

    ### Smoothing for Impair, Imqp at x=d2-d1
    def dIm_minus(x):
        square1 = (x-d21)*(x-d21)
        square2 = (x+d21)*(x+d21)    
        return -x*np.sqrt(d1*d2)*(np.tanh(b*d2)-np.tanh(b*d1))*0.5*np.log((square1+dsm*dsm)*(square2+dsm*dsm)/(square1*square2)) / (4.*d21)

    ### Smoothing for Impair, -Imqp    
    def dIm(x):
        return x*0.5*IP0*((2./np.pi)*np.arctan((1.-x)/dsm) - np.sign(1.-x) + (2./np.pi)*np.arctan((1.+x)/dsm) - np.sign(1.+x))

    ### Smoothing for Impair, Imqp at d2-d1=0
    def dIm_at_0(x):
        x2=x*x
        return -b*x*expb*0.5*np.log((x2+dsm*dsm)/x2)/((expb+1.)*(expb+1.))

    if(symmetric_junction==True):
        def Jpair_correction(x):
            absx = np.maximum(np.abs(x),1.e-5)
            return dRe(absx) + np.sign(x)*1j*(dIm(absx) + dIm_at_0(absx))

        def Jqp_correction(x):
            absx = np.maximum(np.abs(x),1.e-5)
            return dRe(absx) + np.sign(x)*1j*(-dIm(absx) + dIm_at_0(absx))

    else:
        def Jpair_correction(x):
            absx = np.maximum(np.abs(x),1.e-5)
            return dRe(absx) + dRe_minus(absx) + np.sign(x)*1j*(dIm(absx) + dIm_minus(absx))

        def Jqp_correction(x):
            absx = np.maximum(np.abs(x),1.e-5)
            return dRe(absx) - dRe_minus(absx) + np.sign(x)*1j*(-dIm(absx) + dIm_minus(absx))

    def Jpair_smooth(x):     
        return Jpair(x) + Jpair_correction(x)

    def Jqp_smooth(x):
        return Jqp(x) + Jqp_correction(x)

    return (Jpair_smooth, Jqp_smooth)


#-------------------------------------------------------------------------
# Function description: 
#     Fitting by sum of complex exponentials, see Ref. [1,2].
# Input:
#     x -- array of frequency points prepared for fitting
#     Jpair_data -- complex array prepared for fitting of pair TCAs
#     Jqp_data -- complex array prepared for fitting of quasiparticle TCAs
#     maxNterms # number of complex exponentials
#     thr # ratio of absolute and relative tolerances (equal to $\tau_a/\tau_r$ in Ref.[2])
# Output:
#     (p,A,B)
#     Jpair_model
#     Jqp_model
# References:
#    [1] A. A. Odintsov, V. K. Semenov and A. B. Zorin, IEEE Trans. Magn. 23, 763 (1987).
#    [2] D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, Phys. Rev. B 96, 024515 (2017).
#-------------------------------------------------------------------------
def new_fit(x, Jpair_data, Jqp_data, maxNterms, thr):

    assert len(Jpair_data)==len(x) and len(Jqp_data)==len(x)

    p = np.array([-1.+0.j])
    A = np.zeros(1) + 1j*np.zeros(1)
    B = np.zeros(1) + 1j*np.zeros(1)
    xnew=1. # new term to be added at gap frequency

    Nterms = len(p)
    assert len(A)==Nterms and len(B)==Nterms

    start_time = time.time()

    while(Nterms<maxNterms):

        p = np.append(p,-1.+1j*xnew) # new term at frequency xnew
        A = np.append(A,0.+0.j)
        B = np.append(B,0.+0.j)
    
        Nterms += 1
        print('# Nterms = %d with new term at frequency %f. Calculating...' % (Nterms,xnew))

        param_list = leastsq(residual, to_flat((p,A,B)), args=(x, Jpair_data, Jqp_data, thr), ftol=1.e-3)[0] # least square method

        (p,A,B) = from_flat(param_list)

        Jpair_model = modelJpair((p,A,B), x)
        Jqp_model = modelJqp((p,A,B), x)
        ReJpair_diff = Drel(Jpair_model.real, Jpair_data.real, thr)
        ImJpair_diff = Drel(Jpair_model.imag, Jpair_data.imag, thr)
        ReJqp_diff = Drel(Jqp_model.real, Jqp_data.real, thr)
        ImJqp_diff = Drel(Jqp_model.imag, Jqp_data.imag, thr)
        ind_max = np.argmax( np.concatenate((ReJpair_diff,ImJpair_diff,ReJqp_diff,ImJqp_diff)) )
        xnew = x[ind_max%len(x)]

        display((p,A,B))
        print()

    finish_time=time.time()

    print('# Timing: %.f' % (time.time()-start_time),' seconds')

    return ((p,A,B), Jpair_model, Jqp_model)


###-----------------------------------------------
### Supplementary routines for new_fit
###-----------------------------------------------

def display(pAB):
    (p,A,B) = pAB
    print("p=",repr(p))
    print("A=",repr(A))
    print("B=",repr(B))  

# Relative difference with threshold thr
def Drel(X,Xref,thr):
    return np.abs(X-Xref)/np.maximum(thr,np.abs(Xref))

### Jpair model
def modelJpair(pAB, w):
    (p,A,B) = pAB
    rep = p.real
    imp = p.imag
    rea = A.real
    ima = A.imag

    sum=0.0;
    for n in range(len(p)):         
        ### Differs from modelJqp by w->(-w)
        sum-=(-rea[n]*abs(rep[n])+ima[n]*imp[n]-1j*w*rea[n])/(rep[n]**2 + imp[n]**2 - w**2 + 2.*1j*w*abs(rep[n]))

    return sum

### Jqp model
def modelJqp(pAB, w):
    (p,A,B) = pAB
    rep = p.real
    imp = p.imag
    reb = B.real
    imb = B.imag

    sum=0.0;
    for n in range(len(p)):
        sum-=(-reb[n]*abs(rep[n])+imb[n]*imp[n]+1j*w*reb[n])/(rep[n]**2 + imp[n]**2 - w**2 - 2.*1j*w*abs(rep[n]))

    return 1j*w + sum

def to_flat(pAB):
    (p,A,B) = pAB
    pflat = np.stack((p.real,p.imag),axis=1).flatten()
    Aflat = np.stack((A.real,A.imag),axis=1).flatten()
    Bflat = np.stack((B.real,B.imag),axis=1).flatten()
    return np.concatenate((pflat,Aflat,Bflat))

def from_flat(param_list):
    Nterms = int(len(param_list)/6)
    p = -np.abs(param_list[0:2*Nterms:2]) + 1j*param_list[1:2*Nterms:2] # ensure Re[p] < 0
    A = param_list[2*Nterms:4*Nterms:2] + 1j*param_list[1+2*Nterms:4*Nterms:2]
    B = param_list[4*Nterms:6*Nterms:2] + 1j*param_list[1+4*Nterms:6*Nterms:2]
    return (p,A,B)

# Residual for least square optimization
def residual(param_list, x, Jpair_data, Jqp_data, thr):

    (p,A,B) = from_flat(param_list)
    Jpair_model = modelJpair((p,A,B), x)
    Jqp_model = modelJqp((p,A,B), x)

    ReJpair_diff = Drel(Jpair_model.real, Jpair_data.real, thr)
    ImJpair_diff = Drel(Jpair_model.imag, Jpair_data.imag, thr)
    ReJqp_diff = Drel(Jqp_model.real, Jqp_data.real, thr)
    ImJqp_diff = Drel(Jqp_model.imag, Jqp_data.imag, thr)

    return np.concatenate((ReJpair_diff,ImJpair_diff,ReJqp_diff,ImJqp_diff))

#-----------------------------------------------------------------------------------------------
# Function description:
#     Save fit parameters p,A,B to file 
# Input:
#     pAB -- tuple of parameter arrays p,A,B
#     filename where parameters are to be saved
#-----------------------------------------------------------------------------------------------
def parsave(pAB,filename):
    (p,A,B) = pAB 
    np.savetxt(filename,np.c_[p.real,p.imag,A.real,A.imag,B.real,B.imag],fmt='%10.6f')
    print('\n# Parameters saved to file '+filename)
    return

#-----------------------------------------------------------------------------------------------
# Function description:
#     Return fitted TCAs for given parameters 
# Input:
#     pAB -- tuple of parameter arrays p,A,B
# Returns:
#     Tuple of fit functions for pair and quasiparticle TCAs
#-----------------------------------------------------------------------------------------------
def tca_fit(pAB):
    def Jpair_fit(w):
        return modelJpair(pAB, w)
    def Jqp_fit(w):
        return modelJqp(pAB, w)
    return (Jpair_fit, Jqp_fit)

#-----------------------------------------------------------------------------------------------
# Function description:
#     Load TCA fit parameters from file 
# Input:
#     filename with fit parameters p,A,B 
# Returns:
#     Tuple of fit parameter arrays p,A,B
#-----------------------------------------------------------------------------------------------
def load_fit_parameters(filename):
    pars = np.loadtxt(filename)
    Nterms = pars.shape[0]
    p = pars[:,0] + 1j*pars[:,1]
    A = pars[:,2] + 1j*pars[:,3]
    B = pars[:,4] + 1j*pars[:,5]
    return (p,A,B)

#-----------------------------------------------------------------------------------------------
# Function description:
#     Load TCA fit parameters from file and calculate data points at x
# Input:
#     x -- array of frequency points 
#     filename with fit parameters p,A,B 
# Returns:
#     Tuple of fit parameter arrays p,A,B
#     Data points for pair TCA at values of x
#     Data points for quasiparticle TCA at values of x
#-----------------------------------------------------------------------------------------------
def load_fit(x, filename):
    pars = np.loadtxt(filename)
    Nterms = pars.shape[0]
    p = pars[:,0] + 1j*pars[:,1]
    A = pars[:,2] + 1j*pars[:,3]
    B = pars[:,4] + 1j*pars[:,5]
    Jpair_x = modelJpair((p,A,B), x)
    Jqp_x = modelJqp((p,A,B), x)
    return ((p,A,B), Jpair_x, Jqp_x)
